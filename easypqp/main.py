import ast
from email.policy import default
import os
import pkg_resources
import click
import sqlite3
import pandas as pd
from shutil import copyfile
from .convert import conversion, basename_spectralfile
from .library import generate
from .unimoddb import unimod_filter
from easypqp import pkg_unimod_db
from .targetedfileconverter import TargetedFileConverter

try:
    from pyprophet.data_handling import transform_pi0_lambda
except ModuleNotFoundError:
    transform_pi0_lambda = None

@click.group(chain=True)
@click.version_option()
def cli():
    """
    EasyPQP: Simple library generation for OpenSWATH using MSFragger & PyProphet

    Visit https://www.openswath.org for usage instructions and help.
    """


# https://stackoverflow.com/a/47730333
class PythonLiteralOption(click.Option):
    def type_cast_value(self, ctx, value):
        if not isinstance(value, str):  # required for Click>=8.0.0
            return value
        try:
            return ast.literal_eval(value)
        except Exception:
            raise click.BadParameter(value)


# EasyPQP Convert
@cli.command()
@click.option('--pepxml', 'pepxmlfile', required=True, type=click.Path(exists=True), help='The input MSFragger TSV file.')
@click.option('--spectra', 'spectralfile', required=True, type=click.Path(exists=True), help='The input mzXML or MGF (timsTOF only) file.')
@click.option('--unimod', 'unimodfile', required=False, type=click.Path(exists=True), help='The input UniMod XML file.')
@click.option('--psms', 'psmsfile', required=False, type=click.Path(exists=False), help='Output PSMs file.')
@click.option('--peaks', 'peaksfile', required=False, type=click.Path(exists=False), help='Output peaks file.')
@click.option('--exclude-range', 'exclude_range_str', default="-1.5,3.5", show_default=True, required=False, type=str, help='massdiff in this range will not be mapped to UniMod.')
@click.option('--max_delta_unimod', default=0.02, show_default=True, type=float, help='Maximum delta mass (Dalton) for UniMod annotation.')
@click.option('--max_delta_ppm', default=15, show_default=True, type=float, help='Maximum delta mass (PPM) for annotation.')
@click.option('--enable_unannotated/--no-enable_unannotated', default=False, show_default=True, help='Enable mapping uf unannotated delta masses.')
@click.option('--enable_massdiff/--no-enable_massdiff', default=False, show_default=True, help='Enable mapping uf mass differences reported by legacy search engines.')
@click.option('--fragment_types', default="['b','y']", show_default=True, cls=PythonLiteralOption, help='Allowed fragment ion types (a,b,c,x,y,z).')
@click.option('--fragment_charges', default='[1,2,3,4]', show_default=True, cls=PythonLiteralOption, help='Allowed fragment ion charges.')
@click.option('--enable_specific_losses/--no-enable_specific_losses', default=False, show_default=True, help='Enable specific fragment ion losses.')
@click.option('--enable_unspecific_losses/--no-enable_unspecific_losses', default=False, show_default=True, help='Enable unspecific fragment ion losses.')
@click.option('--subsample_fraction', default=1.0, show_default=True, type=float, help='Data fraction used for subsampling.')
@click.option('--max_psm_pep', default=0.5, show_default=True, type=float, help='Maximum posterior error probability (PEP) for a PSM')
def convert(pepxmlfile, spectralfile, unimodfile, psmsfile, peaksfile, exclude_range_str, max_delta_unimod, max_delta_ppm, enable_unannotated, enable_massdiff, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, subsample_fraction, max_psm_pep):
    """
    Convert pepXML files for EasyPQP
    """

    if unimodfile is None:
        unimodfile = pkg_resources.resource_filename('easypqp', 'data/unimod.xml')

    run_id = basename_spectralfile(spectralfile)
    if psmsfile is None:
        psmsfile = run_id + ".psmpkl"
    if peaksfile is None:
        peaksfile = run_id + ".peakpkl"

    temp = exclude_range_str.split(',')
    exclude_range = [float(temp[0]), float(temp[1])]

    click.echo("Info: Converting %s." % pepxmlfile)
    psms, peaks = conversion(pepxmlfile, spectralfile, unimodfile, exclude_range, max_delta_unimod, max_delta_ppm, enable_unannotated, enable_massdiff, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, max_psm_pep)

    psms.to_pickle(psmsfile)
    click.echo("Info: PSMs successfully converted and stored in %s." % psmsfile)

    peaks.to_pickle(peaksfile)
    click.echo("Info: Peaks successfully converted and stored in %s." % peaksfile)

# EasyPQP Library
@cli.command()
@click.argument('infiles', nargs=-1, type=click.Path(exists=True))
@click.option('--out', 'outfile', required=True, type=click.Path(exists=False), help='Output TSV peptide query parameter file.')
@click.option('--psmtsv', 'psmtsv', required=False, type=click.Path(exists=False), help='psm.tsv file from Philosopher.')
@click.option('--peptidetsv', 'peptidetsv', required=False, type=click.Path(exists=False), help='peptide.tsv file from Philosopher.')
@click.option('--perform_rt_calibration', required=False, type=bool, help='Whether to perform RT calibration', default=True, show_default=True)
@click.option('--rt_reference', 'rt_referencefile', required=False, type=click.Path(exists=True), help='Optional iRT/CiRT reference file.')
@click.option('--rt_reference_run_path', 'rt_reference_run_path', default='easypqp_rt_reference_run.tsv', show_default=True, required=False, type=click.Path(exists=False), help='Writes reference run RT file, if RT reference file is not provided.')
@click.option('--rt_filter', 'rt_filter', required=False, type=str, help='Optional tag to filter candidate RT reference runs.')
@click.option('--perform_im_calibration', required=False, type=bool, help='Whether to perform IM calibration', default=True, show_default=True)
@click.option('--im_reference', 'im_referencefile', required=False, type=click.Path(exists=True), help='Optional IM reference file.')
@click.option('--im_reference_run_path', 'im_reference_run_path', default='easypqp_im_reference_run.tsv', show_default=True, required=False, type=click.Path(exists=False), help='Writes reference run IM file, if IM reference file is not provided.')
@click.option('--im_filter', 'im_filter', required=False, type=str, help='Optional tag to filter candidate IM reference runs.')
@click.option('--psm_fdr_threshold', default=0.01, show_default=True, type=float, help='PSM FDR threshold.')
@click.option('--peptide_fdr_threshold', default=0.01, show_default=True, type=float, help='Peptide FDR threshold.')
@click.option('--protein_fdr_threshold', default=0.01, show_default=True, type=float, help='Protein FDR threshold.')
@click.option('--rt_lowess_fraction', default=0.05, show_default=True, type=float, help='Fraction of data points to use for RT lowess regression. If set to 0, cross validation is used.')
@click.option('--rt_psm_fdr_threshold', default=0.001, show_default=True, type=float, help='PSM FDR threshold used for RT alignment.')
@click.option('--im_lowess_fraction', default=0.05, show_default=True, type=float, help='Fraction of data points to use for IM lowess regression. If set to 0, cross validation is used.')
@click.option('--im_psm_fdr_threshold', default=0.001, show_default=True, type=float, help='PSM FDR threshold used for IM alignment.')
@click.option('--pi0_lambda', default=[0.1,0.5,0.05], show_default=True, type=(float, float, float), help='Use non-parametric estimation of p-values. Either use <START END STEPS>, e.g. 0.1, 1.0, 0.1 or set to fixed value, e.g. 0.4, 0, 0.', callback=transform_pi0_lambda)
@click.option('--peptide_plot', 'peptide_plot_path', default="easypqp_peptide_report.pdf", show_default=True, required=True, type=click.Path(exists=False), help='Output peptide-level PDF report.')
@click.option('--protein_plot', 'protein_plot_path', default="easypqp_protein_report.pdf", show_default=True, required=True, type=click.Path(exists=False), help='Output protein-level PDF report.')
@click.option('--min_peptides', default=5, show_default=True, type=int, help='Minimum peptides required for successful alignment.')
@click.option('--proteotypic/--no-proteotypic', show_default=True, default=True, help='Use only proteotypic, unique, non-shared peptides.')
@click.option('--consensus/--no-consensus', show_default=True, default=True, help='Generate consensus instead of best replicate spectra.')
@click.option('--nofdr/--no-fdr-filtering', show_default=True, default=False, help='Do not reassess or filter by FDR, as library was already provided using customized FDR filtering.')
def library(infiles, outfile, psmtsv, peptidetsv, perform_rt_calibration, rt_referencefile, rt_reference_run_path, rt_filter, perform_im_calibration, im_referencefile, im_reference_run_path, im_filter, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, rt_lowess_fraction, rt_psm_fdr_threshold, im_lowess_fraction, im_psm_fdr_threshold, pi0_lambda, peptide_plot_path, protein_plot_path, min_peptides, proteotypic, consensus, nofdr):
    """
    Generate EasyPQP library
    """

    generate(infiles, outfile, psmtsv, peptidetsv, perform_rt_calibration, rt_referencefile, rt_reference_run_path, rt_filter, perform_im_calibration, im_referencefile, im_reference_run_path, im_filter, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, rt_lowess_fraction, rt_psm_fdr_threshold, im_lowess_fraction, im_psm_fdr_threshold, pi0_lambda, peptide_plot_path, protein_plot_path, min_peptides, proteotypic, consensus, nofdr)
    click.echo("Info: Library successfully generated.")

# EasyPQP Reduce
@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Input PQP file.')
@click.option('--out', 'outfile', required=True, type=click.Path(exists=False), help='Output PQP file.')

@click.option('--bins', default=10, show_default=True, type=int, help='Number of bins to fill along gradient.')
@click.option('--peptides', default=5, show_default=True, type=int, help='Number of peptides to sample.')
def reduce(infile, outfile, bins, peptides):
    """
    Reduce PQP files for OpenSWATH linear and non-linear alignment
    """

    # Define outfile
    if outfile is None:
        outfile = infile
    else:
        copyfile(infile, outfile)
        outfile = outfile

    con = sqlite3.connect(outfile)

    anchor_candidates = pd.read_sql('SELECT PRECURSOR.ID AS PRECURSOR_ID, LIBRARY_RT FROM PRECURSOR WHERE PRECURSOR.DECOY == 0;', con)

    anchor_candidates['BIN'] = pd.cut(anchor_candidates['LIBRARY_RT'], bins=bins, right=False, labels=False)

    anchors = anchor_candidates.groupby('BIN').head(peptides).reset_index()

    anchors[['PRECURSOR_ID','LIBRARY_RT']].to_sql('temp_anchors', con, index=False)

    # Delete precursors
    con.execute('DELETE FROM PRECURSOR WHERE ID NOT IN (SELECT PRECURSOR_ID FROM temp_anchors)')

    # Delete transitions
    con.execute('DELETE FROM TRANSITION_PRECURSOR_MAPPING WHERE PRECURSOR_ID NOT IN (SELECT PRECURSOR_ID FROM temp_anchors)')
    con.execute('DELETE FROM TRANSITION WHERE ID NOT IN (SELECT TRANSITION_ID FROM TRANSITION_PRECURSOR_MAPPING)')

    # Delete peptides and proteins
    con.execute('DELETE FROM PRECURSOR_PEPTIDE_MAPPING WHERE PRECURSOR_ID NOT IN (SELECT PRECURSOR_ID FROM temp_anchors)')
    con.execute('DELETE FROM PEPTIDE WHERE ID NOT IN (SELECT PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING)')
    con.execute('DELETE FROM PEPTIDE_PROTEIN_MAPPING WHERE PEPTIDE_ID NOT IN (SELECT PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING)')
    con.execute('DELETE FROM PROTEIN WHERE ID NOT IN (SELECT PROTEIN_ID FROM PEPTIDE_PROTEIN_MAPPING)')

    # Delete tables
    con.execute('DROP TABLE temp_anchors;')

    con.commit()

    # Clean file
    con.execute('VACUUM;')

    # Close connection to file
    con.close()

    click.echo("Info: Library successfully processed and stored in %s." % outfile)

# Parameter transformation functions
def transform_comma_string_to_list(ctx, param, value):
    if value is not None:
        str_list = value.split(",")
        return str_list
    else:
        return None 

# EasyPQP UniMod Database Filtering
@cli.command()
@click.option('--in', 'infile', required=False, default=pkg_unimod_db, show_default=True, type=click.Path(exists=True), help='Input UniMod XML file.')
@click.option('--out', 'outfile', required=False, default="unimod_ipf.xml", show_default=True, type=click.Path(exists=False), help='Output Filtered UniMod XML file.')
@click.option('--ids', 'accession_ids', default='1,2,4,5,7,21,26,27,28,34,35,36,40,121,122,259,267,299,354', show_default=True, type=str, help='UniMod record ids to filter for, i.e. 1,2,4,21.', callback=transform_comma_string_to_list)
@click.option('--sites', 'site_specificity', default=None, show_default=True, type=str, help="""Optional further restriction for specificity, i.e. [n,],M,nK[,QN,STY,*,*,*,EDcRK,WM,RK,Y,K,[TKnS,K,R,EK,Y]. Ensure, you match the sites you want to restrict per unimod.\b\n\nFor example, if --ids=1,21,35, then you should have the following for --sites=n,STY,M. This will restrict acetylation for any N-Term, phosphorylation for serine, threonine, and tyrosine, and oxidation for methionine.
\b\n\nValid Sites:\n\n* - wildcard, will not restrict for any specificty for corresponding UniMod entry.\n\n[ - Protein N-Term\n\n] - Protein C-Term\n\nn - Any N-Term\n\nc - Any C-Term\n\nAmino Acid Letter - A valid amino acid one letter code.\n\n""", callback=transform_comma_string_to_list)
def filter_unimod(infile, outfile, accession_ids, site_specificity):
    """
    Reduce UniMod XML Database file
    """
    unimod_filter(infile, outfile, accession_ids, site_specificity)

# EasyPQP TargetedFileConverter
@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Transition list to convert.')
@click.option('--in_type', default=None, show_default=True, type=str, help='Input file type. Default: None, will be determined from input file. Valid formats: ["tsv", "mrm" ,"pqp", "TraML", "parquet"]')
@click.option('--out', 'outfile', required=False, default="library.pqp", show_default=True, type=click.Path(exists=False), help='Output file to be converted to.')
@click.option('--out_type', default=None, show_default=True, type=str, help='Output file type. Default: None, will be determined from output file. Valid formats: ["tsv", "pqp", "TraML"]')
@click.option('--legacy_traml_id/--no-legacy_traml_id', show_default=True, default=True, help='PQP to TraML: Should legacy TraML IDs be used?')
def targeted_file_converter(infile, in_type, outfile, out_type, legacy_traml_id):
    """
    Convert different spectral libraries / transition files for targeted proteomics and metabolomics analysis.

    Can convert multiple formats to and from TraML (standardized transition format). The following formats are supported:\b\n\n

        - @ref OpenMS::TraMLFile "TraML" \b\n
        - @ref OpenMS::TransitionTSVFile "OpenSWATH TSV transition lists" \b\n
        - @ref OpenMS::TransitionPQPFile "OpenSWATH PQP SQLite files" \b\n
        - SpectraST MRM transition lists \b\n
        - Skyline transition lists \b\n
        - Spectronaut transition lists \b\n
        - Parquet transition lists \b\n

    """
    converter = TargetedFileConverter(infile, outfile, in_type, out_type, legacy_traml_id)
    converter.convert()


# EasyPQP OpenSwathAssayGenerator
@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help="Input file (valid formats: 'tsv', 'mrm', 'pqp', 'TraML')")
@click.option('--in_type', default=None, show_default=True, type=str, help='Input file type. Default: None, will be determined from file extension or content. Valid formats: ["tsv", "mrm" ,"pqp", "TraML"]')
@click.option('--out', 'outfile', required=True, type=click.Path(exists=False), help="Output file (valid formats: 'tsv', 'pqp', 'TraML')")
@click.option('--out_type', default=None, show_default=True, type=str, help='Output file type. Default: None, will be determined from file extension or content. Valid formats: ["tsv", "mrm" ,"pqp", "TraML"]')
@click.option('--min_transitions', default=6, required=False, show_default=True, type=int, help='Minimal number of transitions')
@click.option('--max_transitions', default=6, required=False, show_default=True, type=int, help='Maximal number of transitions')
@click.option('--allowed_fragment_types', required=False, default='b,y', show_default=True, type=str, help='Allowed fragment types')
@click.option('--allowed_fragment_charges', required=False, default='1,2,3,4', show_default=True, type=str, help='Allowed fragment charge states')
@click.option('--enable_detection_specific_losses', required=False, type=bool, help='Set this flag if specific neutral losses for detection fragment ions should be allowed', default=False)
@click.option('--enable_detection_unspecific_losses', required=False, type=bool, help='Set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for detection fragment ions should be allowed', default=False)
@click.option('--precursor_mz_threshold', default=0.025, show_default=True, required=False, type=float, help='MZ threshold in Thomson for precursor ion selection')
@click.option('--precursor_lower_mz_limit', default=400.0, show_default=True, required=False, type=float, help='Lower MZ limit for precursor ions')
@click.option('--precursor_upper_mz_limit', default=1200.0, show_default=True, required=False, type=float, help='Upper MZ limit for precursor ions')
@click.option('--product_mz_threshold', default=0.025, show_default=True, required=False, type=float, help='MZ threshold in Thomson for fragment ion annotation')
@click.option('--product_lower_mz_limit', default=350.0, show_default=True, required=False, type=float, help='Lower MZ limit for fragment ions')
@click.option('--product_upper_mz_limit', default=2000.0, show_default=True, required=False, type=float, help='Upper MZ limit for fragment ions')
@click.option('--swath_windows_file', 'swath_win', required=False, type=click.Path(exists=False), help="Tab separated file containing the SWATH windows for exclusion of fragment ions falling into the precursor isolation window: lower_offset upper_offset \newline400 425 \newline ... Note that the first line is a header and will be skipped. (valid formats: 'txt')")
@click.option('--unimod_file', 'unimod', required=False, type=click.Path(exists=False), help="(Modified) Unimod XML file (http://www.unimod.org/xml/unimod.xml) describing residue modifiability (valid formats: 'xml')")
@click.option('--enable_ipf', required=False, type=bool, help="IPF: set this flag if identification transitions should be generated for IPF. Note: Requires setting 'unimod_file", default=False)
@click.option('--max_num_alternative_localizations', default=10000, required=False, show_default=True, type=int, help='IPF: maximum number of site-localization permutations')
@click.option('--disable_identification_ms2_precursors', required=False, type=bool, help="IPF: set this flag if MS2-level precursor ions for identification should not be allowed for extraction of the precursor signal from the fragment ion data (MS2-level).", default=False)
@click.option('--disable_identification_specific_losses', required=False, type=bool, help="IPF: set this flag if specific neutral losses for identification fragment ions should not be allowed", default=False)
@click.option('--enable_identification_unspecific_losses', required=False, type=bool, help="IPF: set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for identification fragment ions should be allowed", default=False)
@click.option('--enable_swath_specifity', required=False, type=bool, help="IPF: set this flag if identification transitions without precursor specificity (i.e. across whole precursor isolation window instead of precursor MZ) should be generated.", default=False)
def openswath_assay_generator(infile, in_type, outfile, out_type, min_transitions, max_transitions, allowed_fragment_types, allowed_fragment_charges, enable_detection_specific_losses, enable_detection_unspecific_losses, precursor_mz_threshold, precursor_lower_mz_limit, precursor_upper_mz_limit, product_mz_threshold, product_lower_mz_limit, product_upper_mz_limit, swath_win, unimod, enable_ipf, max_num_alternative_localizations, disable_identification_ms2_precursors, disable_identification_specific_losses, enable_identification_unspecific_losses, enable_swath_specifity):
    """
    Generates filtered and optimized assays.

    """
    print("hello")
