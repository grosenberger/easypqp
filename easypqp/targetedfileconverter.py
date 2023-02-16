
import os
import pandas as pd
import pyopenms as po
import ctypes
import click

class TargetedFileConverter:
    '''
    TargetedFileConverter

    Converts different spectral libraries / transition files for targeted proteomics and metabolomics analysis.
  
    Can convert multiple formats to and from TraML (standardized transition format). The following formats are supported:
    <ul>
        <li> @ref OpenMS::TraMLFile "TraML" </li>
        <li> @ref OpenMS::TransitionTSVFile "OpenSWATH TSV transition lists" </li>
        <li> @ref OpenMS::TransitionPQPFile "OpenSWATH PQP SQLite files" </li>
        <li> SpectraST MRM transition lists </li>
        <li> Skyline transition lists </li>
        <li> Spectronaut transition lists </li>
        <li> Parquet transition lists </li>
    </ul>
    '''

    def __init__(self, infile, outfile, in_type=None, out_type=None, legacy_traml_id=True):
        self.infile = infile
        self.outfile = outfile

        # Handle types
        if in_type is None:
            in_type = os.path.splitext(self.infile)[1].split('.')[-1]
        self.in_type = in_type
        if out_type is None:
            out_type = os.path.splitext(self.outfile)[1].split('.')[-1]
        self.out_type = out_type

        self.legacy_traml_id = legacy_traml_id
        self.file_types = po.FileTypes()

        self.tr_exp = po.TargetedExperiment()

    def _validate_type(self):
        
        if self.file_types.nameToType(self.in_type) == po.FileType.UNKNOWN and self.in_type!='parquet':
            click.FileError(filename=self.infile, hint=f"Error: Could not determine input file type! {self.infile}")

        if self.file_types.nameToType(self.out_type) == po.FileType.UNKNOWN:
            click.FileError(filename=self.outfile, hint=f"Error: Could not determine output file type! {self.outfile}")
    
    def convert(self):
        
        self._validate_type()
        
        # If input is parquet, need to write out a temporary tsv to consume for conversion
        if self.in_type == 'parquet':
            tr_list = pd.read_parquet(self.infile)
            temp_in_tsv = f"{os.path.splitext(self.infile)[0]}.tsv"
            tr_list.to_csv(temp_in_tsv, sep="\t")
            # Save org infile information
            self.infile_parquet = self.infile
            self.in_type_parquet = self.in_type
            # Overwrite org infile information with TSV information
            self.infile = temp_in_tsv
            self.in_type = "tsv"

        # Read Input into TargetedExperiment
        if self.file_types.nameToType(self.in_type)-1 == po.FileType.TSV or self.file_types.nameToType(self.in_type)-1 == po.FileType.MRM:
            tsv_reader = po.TransitionTSVFile()
            c_in_file = ctypes.create_string_buffer(self.infile.encode())
            tsv_reader.convertTSVToTargetedExperiment(c_in_file.value, self.file_types.nameToType(self.in_type), self.tr_exp)
        
        elif self.file_types.nameToType(self.in_type)-1 == po.FileType.PQP:
            pqp_reader = po.TransitionPQPFile()
            c_in_file = ctypes.create_string_buffer(self.infile.encode())
            pqp_reader.convertPQPToTargetedExperiment(c_in_file.value, self.tr_exp, self.legacy_traml_id)

        elif self.file_types.nameToType(self.in_type)-1 == po.FileType.TRAML:
            traml_reader = po.TraMLFile()
            c_in_file = ctypes.create_string_buffer(self.infile.encode())
            traml_reader.load(c_in_file, self.tr_exp)

        # Write TargetedExperiment to Output
        if self.file_types.nameToType(self.out_type)-1 == po.FileType.TSV:
            tsv_reader = po.TransitionTSVFile()
            c_out_file = ctypes.create_string_buffer(self.outfile.encode())
            self.tr_exp.getPeptides()
            tsv_reader.convertTargetedExperimentToTSV(c_out_file.value, self.tr_exp)

        elif self.file_types.nameToType(self.out_type)-1 == po.FileType.PQP:
            pqp_reader = po.TransitionPQPFile()
            c_out_file = ctypes.create_string_buffer(self.outfile.encode())
            pqp_reader.convertTargetedExperimentToPQP(c_out_file.value, self.tr_exp)

        elif self.file_types.nameToType(self.out_type)-1 == po.FileType.TRAML:
            traml_reader = po.TraMLFile()
            c_out_file = ctypes.create_string_buffer(self.outfile.encode())
            traml_reader.store(c_out_file, self.tr_exp)

        # Clean Up
        if 'in_type_parquet' in dir(self) and self.out_type!='tsv':
            os.remove(self.infile)
            self.infile = self.infile_parquet

        click.echo(f"Info: Finished converting {self.infile} to {self.outfile}")