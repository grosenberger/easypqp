
import os
import pandas as pd
import pyopenms as po
import ctypes
import click
from typing import Union

class TargetedExperiment:
    """
    Class to load and write a OpenMS TargetedExperiment
    """
    def __init__(self, legacy_traml_id: bool=True) -> None:
        self.legacy_traml_id = legacy_traml_id
        self.tr_exp = po.TargetedExperiment()
        self.file_types = po.FileTypes()

    def _validate_type(self, file: str, file_type: str) -> None:
        """Method to ensure filetype is a known OpenMS compatible transition list file type."""
        if self.file_types.nameToType(file_type) == po.FileType.UNKNOWN and file_type!='parquet':
            click.FileError(filename=file, hint=f"Error: Could not determine file type! {file}")
    
    def _get_file_type(self, infile) -> str:
        """Method to get file type extension from file."""
        return os.path.splitext(infile)[-1].split('.')[-1]

    def _get_file_type_id(self, file_type: str) -> int:
        """Method to get file type id as annotated in OpenMS filetype database."""
        return self.file_types.nameToType(file_type)-1
    
    def load_library(self, infile: str, in_type: Union[str, None]=None) -> None:
        """
        Method to load data from input transition list into an OpenMS TargetedExperiment Object

        Parameters:
            infile: (str) input transition list file to load
            in_type: (str|None) input file type. Default: None. Will be inferred from infile
        """
        if in_type is None:
            in_type = self._get_file_type(infile)
        self._validate_type(infile, in_type)
        # Convert infile str to ctype c char
        c_in_file = ctypes.create_string_buffer(infile.encode())
        if self._get_file_type_id(in_type) == po.FileType.TSV or self._get_file_type_id(in_type) == po.FileType.MRM:
            click.echo("Info: Reading TSV transition list data...")
            tsv_reader = po.TransitionTSVFile()
            tsv_reader.convertTSVToTargetedExperiment(c_in_file.value, self._get_file_type_id(in_type), self.tr_exp)
            tsv_reader.validateTargetedExperiment(self.tr_exp)
        
        elif self._get_file_type_id(self.in_type) == po.FileType.PQP:
            click.echo("Info: Reading PQP transition list data...")
            pqp_reader = po.TransitionPQPFile()
            pqp_reader.convertPQPToTargetedExperiment(c_in_file.value, self.tr_exp, self.legacy_traml_id)
            pqp_reader.validateTargetedExperiment(self.tr_exp)

        elif self._get_file_type_id(in_type) == po.FileType.TRAML:
            click.echo("Info: Reading TraML transition list data...")
            traml_reader = po.TraMLFile()
            traml_reader.load(c_in_file.value, self.tr_exp)
        
        click.echo(f"Info: Loaded {len(self.tr_exp.getCompounds())} Compounds, {len(self.tr_exp.getProteins()) } Proteins, {len(self.tr_exp.getPeptides())} Peptides, and {len(self.tr_exp.getTransitions())} Transitions")
    
    def write_library(self, outfile: str, out_type: Union[str, None]=None) -> None:
        """
        Method to write data from an OpenMS TargetedExperiment Object to disk

        Parameters:
            outfile: (str) output transition list file to load
            out_type: (str|None) output file type Default: None. Will be inferred from infile
        """
        if out_type is None:
            out_type = self._get_file_type(outfile)
        self._validate_type(outfile, out_type)
        # Convert outfile str to ctype c char
        c_out_file = ctypes.create_string_buffer(outfile.encode())
        if self._get_file_type_id(out_type) == po.FileType.TSV:
            click.echo("Info: Writing TSV transition list data to disk...")
            tsv_reader = po.TransitionTSVFile()
            self.tr_exp.getPeptides()
            tsv_reader.convertTargetedExperimentToTSV(c_out_file.value, self.tr_exp)

        elif self._get_file_type_id(out_type) == po.FileType.PQP:
            click.echo("Info: Writing PQP transition list data to disk...")
            pqp_reader = po.TransitionPQPFile()
            pqp_reader.convertTargetedExperimentToPQP(c_out_file.value, self.tr_exp)

        elif self._get_file_type_id(out_type) == po.FileType.TRAML:
            click.echo("Info: Writing TraML transition list data to disk...")
            traml_reader = po.TraMLFile()
            traml_reader.store(c_out_file.value, self.tr_exp)



class TargetedFileConverter(TargetedExperiment):
    '''
    TargetedFileConverter

    Converts different spectral libraries / transition files for targeted proteomics and metabolomics analysis.
  
    Can convert multiple formats to and from TraML (standardized transition format). The following formats are supported:
    
        - @ref OpenMS::TraMLFile "TraML" 
        - @ref OpenMS::TransitionTSVFile "OpenSWATH TSV transition lists" 
        - @ref OpenMS::TransitionPQPFile "OpenSWATH PQP SQLite files" 
        - SpectraST MRM transition lists 
        - Skyline transition lists 
        - Spectronaut transition lists 
        - Parquet transition lists 
    '''

    def __init__(self, infile: str, outfile: str="library.pqp", in_type: Union[str, None]=None, out_type: Union[str, None]=None, legacy_traml_id: bool=True) -> None:
        super().__init__(legacy_traml_id)
        self.infile = infile
        self.outfile = outfile

        # Handle types
        if in_type is None:
            in_type = self._get_file_type(self.infile)
        self.in_type = in_type
        if out_type is None:
            out_type = self._get_file_type(self.outfile)
        self.out_type = out_type

    def convert(self) -> None:
        """Method for converting between spectral library formats"""
        # If input is parquet, need to write out a temporary tsv to consume for conversion
        if self.in_type == 'parquet':
            tr_list = pd.read_parquet(self.infile)
            # Write out a temp tsv file for loading into a TargetedExperiment Object
            temp_in_tsv = f"{os.path.splitext(self.infile)[0]}.tsv"
            tr_list.to_csv(temp_in_tsv, sep="\t")
            # Save org infile information
            self.infile_parquet = self.infile
            self.in_type_parquet = self.in_type
            # Overwrite org infile information with TSV information
            self.infile = temp_in_tsv
            self.in_type = "tsv"    

        # Read Input into TargetedExperiment
        self.load_library(self.infile, self.in_type)

        # Write TargetedExperiment to Output
        self.write_library(self.outfile, self.out_type)

        # Clean Up
        if 'in_type_parquet' in dir(self) and self.out_type!='tsv':
            os.remove(self.infile)
            self.infile = self.infile_parquet

        click.echo(f"Info: Finished converting {self.infile} to {self.outfile}")