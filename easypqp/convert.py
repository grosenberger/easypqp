import numpy as np
import pandas as pd
import os

# Unimod parsing
import xml.etree.ElementTree as ET

# mzXML parsing
import pyopenms as po

class unimod:
	def __init__(self, unimod_file):
		self.unimod_file = unimod_file
		self.ptms = self.parse_unimod()

	def parse_unimod(self):
		namespaces = {'umod': "http://www.unimod.org/xmlns/schema/unimod_2"}
		ET.register_namespace('', "http://www.unimod.org/xmlns/schema/unimod_2")
		tree = ET.parse(self.unimod_file)
		root = tree.getroot()

		ptms = {}
		sites = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V','N-term','C-term']
		positions = ['Anywhere','Any N-term','Any C-term','Protein N-term','Protein C-term']

		for site in sites:
			ptms[site] = {}

		for position in positions:
			for site in sites:
				ptms[site][position] = {}

		for modifications in root.findall('.//umod:modifications', namespaces):
			for modification in modifications.findall('.//umod:mod', namespaces):
				for specificity in modification.findall('.//umod:specificity', namespaces):
					ptms[specificity.attrib['site']][specificity.attrib['position']][int(modification.attrib['record_id'])] = float(modification.findall('.//umod:delta', namespaces)[0].attrib['mono_mass'])

		return ptms

	def get_id(self, site, position, delta_mass):
		min_delta = 0.02
		min_id = -1
		for key, value in self.ptms[site][position].items():
			if abs(value - delta_mass) < min_delta:
				min_id = key

		return(min_id)

def match_modifications(peptide):
	modified_peptide = peptide['peptide_sequence']
	modifications = {}
	if "M|" in peptide['modifications']:
		for modification in peptide['modifications'].split('|')[1:]:
			site, mass = modification.split('$')
			modifications[int(site)] = float(mass)

		for site in sorted(modifications, reverse=True):
			record_id = um.get_id(peptide['peptide_sequence'][site], 'Anywhere', modifications[site])

			if record_id == -1:
				raise click.ClickException("Error: Could not annotate site %s (%s) from peptide %s with delta mass %s." % (site+1, peptide['peptide_sequence'][site], peptide['peptide_sequence'], modifications[site]))

			modified_peptide = modified_peptide[:site+1] + "(UniMod:" + str(record_id) + ")" + modified_peptide[site+1:]

	return modified_peptide

def read_mzxml(mzxml_path, scan_ids):
	fh = po.MzXMLFile()
	fh.setLogType(po.LogType.CMD)
	input_map = po.MSExperiment()
	fh.load(mzxml_path, input_map)

	peaks_list = []
	for scan_id in scan_ids:

		spectrum = input_map.getSpectrum(scan_id - 1)

		product_mzs = []
		intensities = []
		for peak in spectrum:
			product_mzs.append(peak.getMZ())
			intensities.append(peak.getIntensity())

		peaks = pd.DataFrame({'product_mz': product_mzs, 'intensity': intensities})
		peaks['precursor_mz'] = spectrum.getPrecursors()[0].getMZ()
		peaks['scan_id'] = scan_id
		peaks_list.append(peaks)

	if len(peaks_list) > 0:
		transitions = pd.concat(peaks_list)
	else:
		transitions = pd.DataFrame({'product_mz': [], 'precursor_mz': [], 'intensity': [], 'scan_id': [], })
	return(transitions)

def conversion(fraggerfile, mzxmlfile, peptidefile, unimodfile):
	fragger_names = ['scan_id','precursor_neutral_mass','retention_time','precursor_charge','rank','peptide_sequence','upstream_aa','downstream_aa','protein_id','matched_fragments','total_matched_fragments','peptide_neutral_mass','mass_difference','number_tryptic_terminii','number_missed_cleavages','modifications','hyperscore','nextscore','intercept_em','slope_em']
	df = pd.read_table(fraggerfile, header=None, names=fragger_names, index_col=False)

	# Initialize UniMod
	um = unimod(unimodfile)

	# Generate UniMod peptide sequence
	click.echo("Info: Matching modifications to UniMod.")
	df['modified_peptide'] = df[['peptide_sequence','modifications']].apply(match_modifications, axis=1)

	# Update protein identifiers and metadata
	click.echo("Info: Matching peptides to proteins.")
	df = df.drop(columns = 'protein_id')
	peptide_index = pd.read_pickle(peptidefile)
	peptides_1 = df.shape[0]
	df = pd.merge(df, peptide_index, on='peptide_sequence')
	peptides_2 = df.shape[0]

	if peptides_1 != peptides_2:
		raise click.ClickException("Error: Peptides from PSMs (%s) don't match peptides after matching with FASTA (%s). Check digestion parameters." % (peptides_1, peptides_2))

	# Append PyProphet columns
	run_id = os.path.splitext(os.path.basename(fraggerfile))[0]
	df['run_id'] = run_id
	df['group_id'] = df['run_id'] + "_" + df['scan_id'].astype(str)
	df['expect'] = np.power(10,(df['intercept_em'] + (df['hyperscore'] * df['slope_em'])))
	df['var_expectscore'] = 0.0 - np.log(df['expect'])
	df['var_deltascore'] = 1.0 - (df['nextscore'] / df['hyperscore'])
	df['var_lengthscore'] = np.sqrt(df['peptide_sequence'].str.len())
	df['var_charge'] = df['precursor_charge']

	# DIA-Umpire quality tiers
	if run_id.endswith("_Q1"):
		df['var_quality'] = 1
	elif run_id.endswith("_Q2"):
		df['var_quality'] = 2
	elif run_id.endswith("_Q3"):
		df['var_quality'] = 3
	else: # DDA data
		df['var_quality'] = 0

	df = df.rename(index=str, columns={'hyperscore': 'main_var_hyperscore'})

	# Generate spectrum dataframe
	print("Info: Processing spectra from file %s." % sys.argv[5])
	peaks = read_mzxml(sys.argv[5], df['scan_id'].unique().tolist())

	return df, peaks
	