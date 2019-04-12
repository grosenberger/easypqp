import click
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
    click.echo("Info: Peptide index successfully converted and stored in %s." % psmsfile)


# EasyPQP Convert
@cli.command()
@click.option('--fragger', 'fraggerfile', required=True, type=click.Path(exists=True), help='The input MSFragger TSV file.')
@click.option('--mzxml', 'mzxmlfile', required=True, type=click.Path(exists=True), help='The input mzXML file.')
@click.option('--unimodxml', 'unimodxmlfile', required=True, type=click.Path(exists=True), help='The input UniMod XML file.')
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

