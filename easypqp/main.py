import click
import sqlite3
import pandas as pd
from shutil import copyfile
from .index import parse_fasta
from .convert import conversion
from .library import generate

@click.group(chain=True)
@click.version_option()
def cli():
    """
    EasyPQP: Simple library generation for OpenSWATH using MSFragger & PyProphet

    Visit https://www.openswath.org for usage instructions and help.
    """

# EasyPQP Index
@cli.command()
@click.option('--fasta', 'fastafile', required=True, type=click.Path(exists=True), help='The input FASTA file.')
@click.option('--pepidx', 'pepidxfile', required=True, type=click.Path(exists=False), help='The input peptide index file.')
def index(fastafile, pepidxfile):
    """
    Generate peptide index for EasyPQP
    """

    click.echo("Info: Indexing %s." % fastafile)
    pepidx = parse_fasta(fastafile)

    pepidx.to_pickle(pepidxfile)
    click.echo("Info: Peptide index successfully converted and stored in %s." % pepidxfile)


# EasyPQP Convert
@cli.command()
@click.option('--fragger', 'fraggerfile', required=True, type=click.Path(exists=True), help='The input MSFragger TSV file.')
@click.option('--mzxml', 'mzxmlfile', required=True, type=click.Path(exists=True), help='The input mzXML file.')
@click.option('--unimod', 'unimodfile', required=True, type=click.Path(exists=True), help='The input UniMod XML file.')
@click.option('--pepidx', 'pepidxfile', required=True, type=click.Path(exists=True), help='The input peptide index file.')
@click.option('--psms', 'psmsfile', required=True, type=click.Path(exists=False), help='Output PSMs file.')
@click.option('--subpsms', 'subpsmsfile', required=True, type=click.Path(exists=False), help='Output subsampled PSMs file.')
@click.option('--peaks', 'peaksfile', required=True, type=click.Path(exists=False), help='Output peaks file.')
@click.option('--subsample_fraction', default=1.0, show_default=True, type=float, help='Data fraction used for subsampling.')
def convert(fraggerfile, mzxmlfile, unimodfile, pepidxfile, psmsfile, subpsmsfile, peaksfile, subsample_fraction):
    """
    Convert MSFragger TSV files for EasyPQP
    """

    click.echo("Info: Converting %s." % fraggerfile)
    psms, peaks = conversion(fraggerfile, mzxmlfile, unimodfile, pepidxfile)
    subpsms = psms.sample(frac=subsample_fraction)

    psms.to_csv(psmsfile, sep="\t", index=False)
    click.echo("Info: PSMs successfully converted and stored in %s." % psmsfile)

    subpsms.to_csv(subpsmsfile, sep="\t", index=False)
    click.echo("Info: Subsampled PSMs successfully converted and stored in %s." % subpsmsfile)

    peaks.to_pickle(peaksfile)
    click.echo("Info: Peaks successfully converted and stored in %s." % peaksfile)

# EasyPQP Library
@cli.command()
@click.argument('infiles', nargs=-1, type=click.Path(exists=True))
@click.option('--psm_fdr_threshold', default=0.01, show_default=True, type=float, help='PSM FDR threshold.')
@click.option('--peptide_fdr_threshold', default=0.01, show_default=True, type=float, help='Peptide FDR threshold.')
@click.option('--protein_fdr_threshold', default=0.01, show_default=True, type=float, help='Protein FDR threshold.')
@click.option('--peptide_plot', 'peptide_plot_path', required=True, type=click.Path(exists=False), help='Output peptide-level PDF report.')
@click.option('--protein_plot', 'protein_plot_path', required=True, type=click.Path(exists=False), help='Output protein-level PDF report.')
def library(infiles, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, peptide_plot_path, protein_plot_path):
    """
    Generate EasyPQP library
    """

    generate(infiles, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, peptide_plot_path, protein_plot_path)
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


