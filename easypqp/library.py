# plotting
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

import click
import os
import re
import operator
import numpy as np
import pandas as pd

# alignment
from sklearn import preprocessing
import statsmodels.api as sm
from scipy.interpolate import interp1d

# error rate estimation
from pyprophet.stats import pemp, qvalue, pi0est
from pyprophet.ipf import compute_model_fdr

# plotting
from scipy.stats import gaussian_kde
from numpy import linspace, concatenate
from seaborn import lmplot

def plot(path, title, targets, decoys):
  plt.figure(figsize=(10, 5))
  plt.subplots_adjust(hspace=.5)

  plt.subplot(121)
  plt.title("group score distributions")
  plt.xlabel("score")
  plt.ylabel("# of groups")
  plt.hist(
      [targets, decoys], 20, color=['g', 'r'], label=['target', 'decoy'], histtype='bar')
  plt.legend(loc=2)

  plt.subplot(122)
  tdensity = gaussian_kde(targets)
  tdensity.covariance_factor = lambda: .25
  tdensity._compute_covariance()
  ddensity = gaussian_kde(decoys)
  ddensity.covariance_factor = lambda: .25
  ddensity._compute_covariance()
  xs = linspace(min(concatenate((targets, decoys))), max(
      concatenate((targets, decoys))), 200)
  plt.title("group score densities")
  plt.xlabel("score")
  plt.ylabel("density")
  plt.plot(xs, tdensity(xs), color='g', label='target')
  plt.plot(xs, ddensity(xs), color='r', label='decoy')
  plt.legend(loc=2)

  plt.suptitle(title)
  plt.savefig(path)
  plt.close()

def peptide_fdr(psms, peptide_fdr_threshold, pi0_lambda, plot_path):
  pi0_method = 'bootstrap'
  pi0_smooth_df = 3
  pi0_smooth_log_pi0 = False
  pfdr = False

  peptides = psms.groupby(['modified_peptide','decoy'])['pp'].max().reset_index()
  targets = peptides[~peptides['decoy']].copy()
  decoys = peptides[peptides['decoy']].copy()

  targets['p_value'] = pemp(targets['pp'], decoys['pp'])
  targets['q_value'] = qvalue(targets['p_value'], pi0est(targets['p_value'], pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0)['pi0'], pfdr)

  plot(plot_path, "global peptide scores", targets['pp'], decoys['pp'])
  
  return targets[targets['q_value'] < peptide_fdr_threshold]['modified_peptide']

def protein_fdr(psms, protein_fdr_threshold, pi0_lambda, plot_path):
  pi0_method = 'bootstrap'
  pi0_smooth_df = 3
  pi0_smooth_log_pi0 = False
  pfdr = False

  proteins = psms.groupby(['protein_id','decoy'])['pp'].max().reset_index()
  targets = proteins[~proteins['decoy']].copy()
  decoys = proteins[proteins['decoy']].copy()

  targets['p_value'] = pemp(targets['pp'], decoys['pp'])
  targets['q_value'] = qvalue(targets['p_value'], pi0est(targets['p_value'], pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0)['pi0'], pfdr)

  plot(plot_path, "global protein scores", targets['pp'], decoys['pp'])
  
  return targets[targets['q_value'] < protein_fdr_threshold]['protein_id']

def process_psms(psms, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, pi0_lambda, peptide_plot_path, protein_plot_path, proteotypic):
  # Append columns
  psms['base_name'] = psms['run_id'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

  # Filter proteotypic peptides
  if proteotypic:
    psms = psms[psms['num_tot_proteins'] == 1].copy()
  else:
    raise click.ClickException("Support for non-proteotypic peptides is not yet implemented.")

  # Prepare PeptideProphet / iProphet results
  if 'q_value' not in psms.columns:
    psms['q_value'] = compute_model_fdr(psms['pep'].values)

  # Confident peptides and protein in global context
  peptides = peptide_fdr(psms, peptide_fdr_threshold, pi0_lambda, peptide_plot_path)
  click.echo("Info: %s modified peptides identified (q-value < %s)" % (len(peptides), peptide_fdr_threshold))
  proteins = protein_fdr(psms, protein_fdr_threshold, pi0_lambda, protein_plot_path)
  click.echo("Info: %s proteins identified (q-value < %s)" % (len(proteins), protein_fdr_threshold))

  # Filter peptides and proteins
  psms = psms[psms['modified_peptide'].isin(peptides)]
  psms = psms[psms['protein_id'].isin(proteins)]

  # Filter PSMs
  psms = psms[psms['q_value'] < psm_fdr_threshold]

  # Remove decoys
  psms = psms[~psms['decoy']]

  click.echo("Info: %s redundant PSMs identified (q-value < %s)" % (psms.shape[0], psm_fdr_threshold))

  return psms

def lowess(run, reference_run, min_peptides, base_name, main_path):
  dfm = pd.merge(run, reference_run[['modified_peptide','precursor_charge','irt']], on=['modified_peptide','precursor_charge'])
  click.echo("Info: Peptide overlap between run and reference: %s." % (dfm.shape[0]))
  if dfm.shape[0] <= min_peptides:
    click.echo("Info: Skipping run because not enough peptides could be found for alignment.")
    return pd.DataFrame()

  # Fit lowess model
  lwf = sm.nonparametric.lowess(dfm['irt'], dfm['retention_time'], frac=.66)
  lwf_x = list(zip(*lwf))[0]
  lwf_y = list(zip(*lwf))[1]
  lwi = interp1d(lwf_x, lwf_y, bounds_error=False, fill_value="extrapolate")

  # Apply lowess model
  run['irt'] = lwi(run['retention_time'])

  # Plot regression
  fig = lmplot(x='retention_time', y='irt', data=dfm, lowess=True)
  fig.savefig(os.path.join(main_path, "easypqp_alignment_" + base_name + ".pdf"))
  plt.close()

  return run

def generate(files, referencefile, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, pi0_lambda, peptide_plot_path, protein_plot_path, min_peptides, proteotypic):
  # Parse input arguments
  psm_files = []
  spectra = []

  for file in files:
    if 'psms' in file:
      psm_files.append(file)
    if 'peakpkl' in file:
      spectra.append(file)

  if len(psm_files) == 0:
    raise click.ClickException("No PSMs files present. Need to have tag 'psms' in filename.")

  if len(spectra) == 0:
    raise click.ClickException("No spectrum files present. Need to have tag 'peakpkl' in filename.")

  # Read all PSM files
  psms_list = []
  for psm_file in psm_files:
    click.echo("Info: Reading file %s." % psm_file)
    psms_list.append(pd.read_csv(psm_file, index_col=False, sep='\t'))
  psms = pd.concat(psms_list).reset_index(drop=True)
  psms['pp'] = 1-psms['pep']

  # Process PSMs
  pepid = process_psms(psms, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, pi0_lambda, peptide_plot_path, protein_plot_path, proteotypic)

  # Get main path for figures
  main_path = os.path.dirname(os.path.abspath(peptide_plot_path))

  # Generate set of best replicate identifications per run
  pepidr = pepid.loc[pepid.groupby(['base_name','modified_peptide','precursor_charge'])['pp'].idxmax()].sort_index()

  # Prepare reference iRT list
  if referencefile is not None:
    # Read reference file if present
    reference_run = pd.read_csv(referencefile, index_col=False, sep='\t')
    align_runs = pepidr
    if not set(['modified_peptide','precursor_charge','irt']).issubset(reference_run.columns):
      raise click.ClickException("Reference iRT file has wrong format. Requires columns 'modified_peptide', 'precursor_charge' and 'irt'.")
    if reference_run.shape[0] < 10:
      raise click.ClickException("Reference iRT file has too few data points. Requires at least 10.")
  else:
    # Select reference run
    pepidr_stats = pepidr.groupby('base_name')[['modified_peptide']].count().reset_index()
    click.echo(pepidr_stats)
    reference_run_base_name = pepidr_stats.loc[pepidr_stats['modified_peptide'].idxmax()]['base_name']

    reference_run = pepidr[pepidr['base_name'] == reference_run_base_name].copy()
    align_runs = pepidr[pepidr['base_name'] != reference_run_base_name]

    # Normalize RT of reference run
    min_max_scaler = preprocessing.MinMaxScaler()
    reference_run['irt'] = min_max_scaler.fit_transform(reference_run[['retention_time']])*100

  # Normalize RT of all runs against reference
  aligned_runs = align_runs.groupby('base_name').apply(lambda x: lowess(x, reference_run, min_peptides, x.name, main_path))
  pepida = aligned_runs
  if referencefile is None:
    pepida = pd.concat([reference_run, aligned_runs], sort=True).reset_index(drop=True)

  # Generate set of non-redundant global best replicate identifications
  pepidb = pepida.loc[pepida.groupby(['modified_peptide','precursor_charge'])['pp'].idxmax()].sort_index()

  # Prepare ID mzML pairing
  peak_files = pd.DataFrame({'path': spectra})
  peak_files['base_name'] = peak_files['path'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

  # Parse mzXML to retrieve peaks and store results in peak files
  for idx, peak_file in peak_files.iterrows():
    click.echo("Info: Parsing file %s." % peak_file['path'])
    meta_run = pepida[pepida['base_name'] == peak_file['base_name']]
    meta_global = pepidb[pepidb['base_name'] == peak_file['base_name']]
    peaks = pd.read_pickle(peak_file['path'])
    
    # Generate run-specific PQP files for OpenSWATH alignment
    if "_Q1" in peak_file['base_name']:
      run_pqp = pd.merge(meta_run, peaks, on='scan_id')[['precursor_mz','product_mz','intensity','irt','protein_id','peptide_sequence','modified_peptide','precursor_charge']]
      run_pqp.columns = ['PrecursorMz','ProductMz','LibraryIntensity','NormalizedRetentionTime','ProteinId','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
      run_pqp['PrecursorCharge'] = run_pqp['PrecursorCharge'].astype(int)
      run_pqp_path = os.path.splitext(peak_file['path'])[0]+"_run_peaks.tsv"
      run_pqp.to_csv(run_pqp_path, sep="\t", index=False)

    # Generate global non-redundant PQP files
    global_pqp = pd.merge(meta_global, peaks, on='scan_id')[['precursor_mz','product_mz','intensity','irt','protein_id','peptide_sequence','modified_peptide','precursor_charge']]
    global_pqp.columns = ['PrecursorMz','ProductMz','LibraryIntensity','NormalizedRetentionTime','ProteinId','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
    global_pqp['PrecursorCharge'] = global_pqp['PrecursorCharge'].astype(int)
    global_pqp_path = os.path.splitext(peak_file['path'])[0]+"_global_peaks.tsv"
    global_pqp.to_csv(global_pqp_path, sep="\t", index=False)
