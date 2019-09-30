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
import sklearn.isotonic
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

def process_psms(psms, psmtsv, peptidetsv, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, pi0_lambda, peptide_plot_path, protein_plot_path, proteotypic):
  # Append columns
  psms['base_name'] = psms['run_id'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

  if None not in (psmtsv, peptidetsv):
    # Read psm.tsv and peptide.tsv
    peptidetsv_df = pd.read_csv(peptidetsv, index_col=False, sep='\t', usecols=["Peptide", "Gene", "Protein ID"])
    psmtsv_df = pd.read_csv(psmtsv, index_col=False, sep='\t', usecols=["Spectrum", "Peptide"])

    # Filter out PSMs whose peptides are not in peptide.tsv
    psmtsv_df = psmtsv_df[psmtsv_df["Peptide"].isin(peptidetsv_df["Peptide"])]

    # Generate a group_id column
    temp_df = psmtsv_df["Spectrum"].str.split("\\.", expand=True)
    psmtsv_df["group_id"] = temp_df.iloc[:, 0] + "_" + pd.to_numeric(temp_df.iloc[:, -2]).astype(str)

    # Filter psm dataframe
    psms = psms[psms["group_id"].isin(psmtsv_df["group_id"])]

    # Update gene_id and protein_id
    psms = psms.merge(peptidetsv_df, how="left", left_on="peptide_sequence", right_on="Peptide")
    psms.drop(["gene_id", "protein_id"], inplace=True, axis=1)
    psms.rename(columns={"Gene": "gene_id", "Protein ID": "protein_id"}, inplace=True)
    psms["num_tot_proteins"] = 1
    click.echo("Info: %s redundant PSMs identified after filtering with %s and %s" % (psms.shape[0], psmtsv, peptidetsv))
  else:
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

def lowess(run, reference_run, xcol, ycol, lowess_frac, min_peptides, filename, main_path):
  dfm = pd.merge(run, reference_run[['modified_peptide','precursor_charge',ycol]], on=['modified_peptide','precursor_charge'])
  click.echo("Info: Peptide overlap between run and reference: %s." % (dfm.shape[0]))
  if dfm.shape[0] <= min_peptides:
    click.echo("Info: Skipping run because not enough peptides could be found for alignment.")
    return pd.DataFrame()

  # Fit lowess model
  lwf = sm.nonparametric.lowess(dfm[ycol], dfm[xcol], frac=lowess_frac)
  lwf_x = lwf[:, 0]
  ir = sklearn.isotonic.IsotonicRegression() # make the regression strictly increasing
  lwf_y = ir.fit_transform(lwf_x, lwf[:, 1])
  mask = np.concatenate([[True], np.diff(lwf_y) != 0]) # remove non increasing points
  lwi = interp1d(lwf_x[mask], lwf_y[mask], bounds_error=False, fill_value="extrapolate")

  # Apply lowess model
  run[ycol] = lwi(run[xcol])

  # Plot regression
  plt.plot(dfm[xcol], dfm[ycol], 'o')
  plt.plot(run[xcol], run[ycol])
  plt.xlabel(xcol)
  plt.ylabel(ycol)
  plt.savefig(os.path.join(main_path, filename + ".pdf"))
  plt.close()

  return run

def generate(files, outfile, psmtsv, peptidetsv, rt_referencefile, im_referencefile, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, rt_lowess_frac, im_lowess_frac, pi0_lambda, peptide_plot_path, protein_plot_path, min_peptides, proteotypic, consensus):
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

  if peptidetsv is not None and psmtsv is None:
    raise click.ClickException("There is a peptide.tsv but no psm.tsv.")
  elif peptidetsv is None and psmtsv is not None:
    raise click.ClickException("There is a psm.tsv but no peptide.tsv.")

  if None not in (psmtsv, peptidetsv):
    click.echo("Info: There are psm.tsv and peptide.tsv. Will ignore --psm_fdr_threshold, --peptide_fdr_threshold, --protein_fdr_threshold, --pi0_lambda, ----proteotypic, and --no-proteotypic.")

  # Read all PSM files
  psms_list = []
  for psm_file in psm_files:
    click.echo("Info: Reading file %s." % psm_file)
    psm_tab = pd.read_csv(psm_file, index_col=False, sep='\t')
    if psm_tab.shape[0] > 0:
      psms_list.append()
  psms = pd.concat(psms_list).reset_index(drop=True)
  psms['pp'] = 1-psms['pep']

  # Process PSMs
  pepid = process_psms(psms, psmtsv, peptidetsv, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold, pi0_lambda, peptide_plot_path, protein_plot_path, proteotypic)

  # Get main path for figures
  main_path = os.path.dirname(os.path.abspath(peptide_plot_path))

  # Generate set of best replicate identifications per run
  pepidr = pepid.loc[pepid.groupby(['base_name','modified_peptide','precursor_charge'])['pp'].idxmax()].sort_index()

  # Prepare reference IM list
  enable_im = False
  if im_referencefile is not None:
    enable_im = True
    # Read reference file if present
    im_reference_run = pd.read_csv(im_referencefile, index_col=False, sep='\t')
    align_runs = pepidr
    if not set(['modified_peptide','precursor_charge','im']).issubset(im_referencefile.columns):
      raise click.ClickException("Reference IM file has wrong format. Requires columns 'modified_peptide', 'precursor_charge' and 'im'.")
    if im_reference_run.shape[0] < 10:
      raise click.ClickException("Reference IM file has too few data points. Requires at least 10.")
  elif 'ion_mobility' in pepidr.columns:
    enable_im = True
    # Select reference run
    pepidr_stats = pepidr.groupby('base_name')[['modified_peptide']].count().reset_index()
    click.echo(pepidr_stats)
    im_reference_run_base_name = pepidr_stats.loc[pepidr_stats['modified_peptide'].idxmax()]['base_name']

    im_reference_run = pepidr[pepidr['base_name'] == im_reference_run_base_name].copy()

    # Set IM of reference run
    im_reference_run['im'] = im_reference_run['ion_mobility']

  # Prepare reference iRT list
  if rt_referencefile is not None:
    # Read reference file if present
    rt_reference_run = pd.read_csv(rt_referencefile, index_col=False, sep='\t')
    align_runs = pepidr
    if not set(['modified_peptide','precursor_charge','irt']).issubset(rt_reference_run.columns):
      raise click.ClickException("Reference iRT file has wrong format. Requires columns 'modified_peptide', 'precursor_charge' and 'irt'.")
    if rt_reference_run.shape[0] < 10:
      raise click.ClickException("Reference iRT file has too few data points. Requires at least 10.")
  else:
    # Select reference run
    pepidr_stats = pepidr.groupby('base_name')[['modified_peptide']].count().reset_index()
    click.echo(pepidr_stats)
    rt_reference_run_base_name = pepidr_stats.loc[pepidr_stats['modified_peptide'].idxmax()]['base_name']

    rt_reference_run = pepidr[pepidr['base_name'] == rt_reference_run_base_name].copy()

    # Normalize RT of reference run
    min_max_scaler = preprocessing.MinMaxScaler()
    rt_reference_run['irt'] = min_max_scaler.fit_transform(rt_reference_run[['retention_time']])*100

  # Normalize RT of all runs against reference
  aligned_runs = align_runs.groupby('base_name').apply(lambda x: lowess(x, rt_reference_run, 'retention_time', 'irt', rt_lowess_frac, min_peptides, "easypqp_rt_alignment_" + x.name, main_path))

  # Normalize IM of all runs against reference
  if enable_im:
    aligned_runs = aligned_runs.groupby('base_name').apply(lambda x: lowess(x, im_reference_run, 'ion_mobility', 'im', im_lowess_frac, min_peptides, "easypqp_im_alignment_" + x.name, main_path))
    
  pepida = aligned_runs

  # Remove peptides without valid iRT
  pepida = pepida.loc[np.isfinite(pepida['irt'])]

  # Remove peptides without valid IM
  if enable_im:
    pepida = pepida.loc[np.isfinite(pepida['im'])]
  else:
    pepida.loc[:, 'im'] = np.nan

  # Generate set of non-redundant global best replicate identifications
  pepidb = pepida.loc[pepida.groupby(['modified_peptide','precursor_charge'])['pp'].idxmax()].sort_index()

  # Prepare ID mzML pairing
  peak_files = pd.DataFrame({'path': spectra})
  peak_files['base_name'] = peak_files['path'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

  # Parse mzXML to retrieve peaks and store results in peak files
  replicate_pqp = []
  for idx, peak_file in peak_files.iterrows():
    click.echo("Info: Parsing file %s." % peak_file['path'])
    meta_run = pepida[pepida['base_name'] == peak_file['base_name']]
    meta_global = pepidb[pepidb['base_name'] == peak_file['base_name']]
    peaks = pd.read_pickle(peak_file['path'])
    
    # Generate run-specific PQP files for OpenSWATH alignment
    if consensus or ("_Q1" in peak_file['base_name']):
      run_pqp = pd.merge(meta_run, peaks, on=['modified_peptide','precursor_charge','scan_id'])[['precursor_mz','product_mz','fragment','intensity','irt','im','protein_id','gene_id','peptide_sequence','modified_peptide','precursor_charge']]
      run_pqp.columns = ['PrecursorMz','ProductMz','Annotation','LibraryIntensity','NormalizedRetentionTime','PrecursorIonMobility','ProteinId','GeneName','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
      run_pqp['PrecursorCharge'] = run_pqp['PrecursorCharge'].astype(int)
      run_pqp_path = os.path.splitext(peak_file['path'])[0]+"_run_peaks.tsv"
      run_pqp.to_csv(run_pqp_path, sep="\t", index=False)
      if consensus:
        replicate_pqp.append(run_pqp)

    # Generate global non-redundant PQP files
    if not consensus:
      global_pqp = pd.merge(meta_global, peaks, on=['modified_peptide','precursor_charge','scan_id'])[['precursor_mz','product_mz','fragment','intensity','irt','im','protein_id','gene_id','peptide_sequence','modified_peptide','precursor_charge']]
      global_pqp.columns = ['PrecursorMz','ProductMz','Annotation','LibraryIntensity','NormalizedRetentionTime','PrecursorIonMobility','ProteinId','GeneName','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
      global_pqp['PrecursorCharge'] = global_pqp['PrecursorCharge'].astype(int)
      replicate_pqp.append(global_pqp)

  # Aggregate consensus spectra
  pqp = pd.concat(replicate_pqp)
  if consensus:
    pqp_irt = pqp[['ModifiedPeptideSequence','PrecursorCharge','NormalizedRetentionTime','PrecursorIonMobility']].drop_duplicates().groupby(['ModifiedPeptideSequence','PrecursorCharge'])[['NormalizedRetentionTime','PrecursorIonMobility']].median().reset_index()
    pqp_mass = pqp.groupby(['PrecursorMz','ProductMz','Annotation','ProteinId','GeneName','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge'])['LibraryIntensity'].median().reset_index()
    pqp = pd.merge(pqp_mass,pqp_irt, on=['ModifiedPeptideSequence','PrecursorCharge'])

  # Write output TSV file
  pqp.to_csv(outfile, sep="\t", index=False)
