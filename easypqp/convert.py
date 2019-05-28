import numpy as np
import pandas as pd
import os
from statistics import median_low
import click

# Unimod parsing
import xml.etree.ElementTree as ET

# mzXML parsing
import pyopenms as po

class pepxml:
	def __init__(self, pepxml_file, unimod):
		self.pepxml_file = pepxml_file
		self.psms = self.parse_pepxml()
		self.match_unimod(unimod)

	def get(self):
		return(self.psms)

	def match_unimod(self, unimod):
		def match_modifications(um, peptide):
			monomeric_masses = {"A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694, "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146, "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496, "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203, "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841}
			modified_peptide = peptide['peptide_sequence']

			# parse terminal modifications
			nterm_modification = ""
			if peptide['nterm_modification'] is not "":
				nterm_modification = peptide['nterm_modification'] - 1.0078
			cterm_modification = ""
			if peptide['cterm_modification'] is not "":
				nterm_modification = peptide['cterm_modification']

			# parse closed modifications
			modifications = {}
			if "M|" in peptide['modifications']:
				for modification in peptide['modifications'].split('|')[1:]:
					site, mass = modification.split('$')
					delta_mass = float(mass) - monomeric_masses[peptide['peptide_sequence'][int(site)-1]]
					modifications[int(site)] = delta_mass

			# parse open modifications
			oms_sequence = peptide['peptide_sequence']
			for site in modifications.keys():
				oms_sequence = oms_sequence[:site-1] + "_" + oms_sequence[site:]

			oms_modifications, nterm_modification, cterm_modification = um.get_oms_id(oms_sequence, peptide['massdiff'], nterm_modification, cterm_modification)
			modifications = {**modifications, **oms_modifications}

			for site in sorted(modifications, reverse=True):
				record_id = um.get_id(peptide['peptide_sequence'][site-1], 'Anywhere', modifications[site])

				if record_id == -1:
					raise click.ClickException("Error: Could not annotate site %s (%s) from peptide %s with delta mass %s." % (site, peptide['peptide_sequence'][site-1], peptide['peptide_sequence'], modifications[site]))

				modified_peptide = modified_peptide[:site] + "(UniMod:" + str(record_id) + ")" + modified_peptide[site:]

			if nterm_modification is not "":
				record_id_nterm = um.get_id("N-term", 'Any N-term', nterm_modification)
				if record_id_nterm == -1:
					record_id_nterm = um.get_id("N-term", 'Protein N-term', nterm_modification)

				if record_id_nterm == -1:
					raise click.ClickException("Error: Could not annotate N-terminus from peptide %s with delta mass %s." % (peptide['peptide_sequence'], nterm_modification))
				
				modified_peptide = ".(UniMod:" + str(record_id_nterm) + ")" + modified_peptide

			if cterm_modification is not "":
				record_id_cterm = um.get_id("C-term", 'Any C-term', cterm_modification)
				if record_id_cterm == -1:
					record_id_cterm = um.get_id("C-term", 'Protein C-term', cterm_modification)

				if record_id_cterm == -1:
					raise click.ClickException("Error: Could not annotate C-terminus from peptide %s with delta mass %s." % (peptide['peptide_sequence'], cterm_modification))
				
				modified_peptide = modified_peptide + ".(UniMod:" + str(record_id_cterm) + ")"


			return modified_peptide

		self.psms['modified_peptide'] = self.psms[['peptide_sequence','modifications','nterm_modification','cterm_modification','massdiff']].apply(lambda x: match_modifications(unimod, x), axis=1)

	def parse_pepxml(self):
		peptides = []
		namespaces = {'pepxml_ns': "http://regis-web.systemsbiology.net/pepXML"}
		ET.register_namespace('', "http://regis-web.systemsbiology.net/pepXML")
		tree = ET.parse(self.pepxml_file)
		root = tree.getroot()

		for msms_run_summary in root.findall('.//pepxml_ns:msms_run_summary', namespaces):
			base_name = os.path.basename(msms_run_summary.attrib['base_name'])

			# find decoy prefix
			decoy_prefix = ""
			for search_summary in msms_run_summary.findall('.//pepxml_ns:search_summary', namespaces):
				for parameter in search_summary.findall('.//pepxml_ns:parameter', namespaces):
					if parameter.attrib['name'] == 'decoy_prefix':
						decoy_prefix = parameter.attrib['value']

			# go through all spectrum queries
			for spectrum_query in msms_run_summary.findall('.//pepxml_ns:spectrum_query', namespaces):
				index = spectrum_query.attrib['index']
				start_scan = spectrum_query.attrib['start_scan']
				end_scan = spectrum_query.attrib['end_scan']
				assumed_charge = spectrum_query.attrib['assumed_charge']
				retention_time_sec = spectrum_query.attrib['retention_time_sec']

				for search_result in spectrum_query.findall(".//pepxml_ns:search_result", namespaces):
					for search_hit in search_result.findall(".//pepxml_ns:search_hit", namespaces):
						hit_rank = search_hit.attrib['hit_rank']
						massdiff = search_hit.attrib['massdiff']

						# parse peptide and protein information
						peptide = search_hit.attrib['peptide']
						unprocessed_proteins = [search_hit.attrib['protein']]

						for alternative_protein in search_hit.findall('.//pepxml_ns:alternative_protein', namespaces):
							unprocessed_proteins.append(alternative_protein.attrib['protein'])

						# remove decoy results from mixed target/decoy hits
						has_targets = False
						has_decoys = False
						for prot in unprocessed_proteins:
							if decoy_prefix in prot:
								has_decoys = True
							else:
								has_targets = True

						processed_proteins = []
						for prot in unprocessed_proteins:
							if has_targets and has_decoys:
								if decoy_prefix not in prot:
									processed_proteins.append(prot)
							else:
								processed_proteins.append(prot)
						num_tot_proteins = len(processed_proteins)

						is_decoy = False
						if has_decoys and not has_targets:
							is_decoy = True

						proteins = {}
						for prot in processed_proteins:
							# Remove UniProt prefixes if necessary
							if decoy_prefix + "sp|" in prot:
								proteins[decoy_prefix + prot.split("|")[1]] = ""
							elif "sp|" in prot:
								proteins[prot.split("|")[1]] = prot.split("|")[2]
							else:
								proteins[prot] = prot

						protein = ""
						protein_description = ""
						for key in sorted(proteins):
							if protein == "":
								protein = key
								protein_description = proteins[key]
							else:
								protein = protein + ";" + key
								protein_description = protein_description + ";" + proteins[key]

						# parse PTM information
						modifications = "M"
						nterm_modification = ""
						cterm_modification = ""
						for modification_info in search_hit.findall('.//pepxml_ns:modification_info', namespaces):
							if 'mod_nterm_mass' in modification_info.attrib:
								nterm_modification = float(modification_info.attrib['mod_nterm_mass'])
							if 'mod_cterm_mass' in modification_info.attrib:
								cterm_modification = float(modification_info.attrib['mod_cterm_mass'])
							for mod_aminoacid_mass in modification_info.findall('.//pepxml_ns:mod_aminoacid_mass', namespaces):
								modifications = modifications + "|" + mod_aminoacid_mass.attrib['position'] + "$" + mod_aminoacid_mass.attrib['mass']

						# parse search engine score information
						scores = {}
						for search_score in search_hit.findall('.//pepxml_ns:search_score', namespaces):
							scores["var_" + search_score.attrib['name']] = float(search_score.attrib['value'])

						# parse PeptideProphet or iProphet results if available
						for analysis_result in search_hit.findall('.//pepxml_ns:analysis_result', namespaces):
							if analysis_result.attrib['analysis'] == 'interprophet':
								for interprophet_result in analysis_result.findall('.//pepxml_ns:interprophet_result', namespaces):
									scores["pep"] = 1.0 - float(interprophet_result.attrib['probability'])
							elif analysis_result.attrib['analysis'] == 'peptideprophet':
								for peptideprophet_result in analysis_result.findall('.//pepxml_ns:peptideprophet_result', namespaces):
									scores["pep"] = 1.0 - float(peptideprophet_result.attrib['probability'])

						peptides.append({**{'run_id': base_name, 'scan_id': int(start_scan), 'hit_rank': int(hit_rank), 'massdiff': float(massdiff), 'precursor_charge': int(assumed_charge), 'retention_time': float(retention_time_sec), 'peptide_sequence': peptide, 'modifications': modifications, 'nterm_modification': nterm_modification, 'cterm_modification': cterm_modification, 'protein_id': protein, 'protein_description': protein_description, 'num_tot_proteins': num_tot_proteins, 'decoy': is_decoy}, **scores})

		df = pd.DataFrame(peptides)
		return(df)

class unimod:
	def __init__(self, unimod_file, max_delta):
		self.unimod_file = unimod_file
		self.max_delta = max_delta
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
		candidates = {}
		min_id = -1
		for key, value in self.ptms[site][position].items():
			delta_mod = abs(value - float(delta_mass))
			if delta_mod < self.max_delta:
				if key in candidates.keys():
					if delta_mod < candidates[key]:
						candidates[key] = delta_mod
				else:
					candidates[key] = delta_mod

		if len(candidates) > 0:
			min_id = min(candidates, key=candidates.get)

		return(min_id)

	def get_oms_id(self, sequence, massdiff, nterm_modification, cterm_modification):
		record_ids = {}
		for site, aa in enumerate(sequence):
			if aa != "_":
				record_id_site = self.get_id(aa, 'Anywhere', massdiff)
				if record_id_site != -1:
					record_ids[site+1] = record_id_site

		record_id_nterm = -1
		if nterm_modification == "":
			record_id_nterm = self.get_id("N-term", 'Any N-term', massdiff)
			if record_id_nterm == -1:
				record_id_nterm = self.get_id("N-term", 'Protein N-term', massdiff)

		record_id_cterm = -1
		if cterm_modification == "":
			record_id_cterm = self.get_id("C-term", 'Any C-term', massdiff)
			if record_id_cterm == -1:
				record_id_cterm = self.get_id("C-term", 'Protein C-term', massdiff)

		# prefer residual over N-term over C-term modifications
		aamod = {}
		if len(record_ids) > 0:
			aasite = median_low(list(record_ids.keys()))
			aamod[aasite] = massdiff
		elif record_id_nterm != -1:
			nterm_modification = massdiff
		elif record_id_cterm != -1:
			cterm_modification = massdiff

		return aamod, nterm_modification, cterm_modification

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

def conversion(pepxmlfile, mzxmlfile, unimodfile, main_score, max_delta):
	# Initialize UniMod
	um = unimod(unimodfile, max_delta)

	# Parse pepXML
	px = pepxml(pepxmlfile, um)
	psms = px.get()

	# Generate UniMod peptide sequence
	click.echo("Info: Matching modifications to UniMod.")

	# Append PyProphet columns
	run_id = os.path.splitext(os.path.basename(pepxmlfile))[0]
	psms['group_id'] = psms['run_id'] + "_" + psms['scan_id'].astype(str)

	if 'var_expect' in psms.columns:
		psms = psms.rename(index=str, columns={'var_expect': 'expect'})
		psms['var_expectscore'] = 0.0 - np.log(psms['expect'])

	if 'var_nextscore' in psms.columns and 'var_hyperscore' in psms.columns:
		psms = psms.rename(index=str, columns={'var_nextscore': 'nextscore'})
		psms['var_deltascore'] = 1.0 - (psms['nextscore'] / psms['var_hyperscore'])

	# DIA-Umpire quality tiers
	if run_id.endswith("_Q1"):
		psms['quality'] = 1
	elif run_id.endswith("_Q2"):
		psms['quality'] = 2
	elif run_id.endswith("_Q3"):
		psms['quality'] = 3
	else: # DDA data
		psms['quality'] = 0

	if main_score not in psms.columns:
		raise click.ClickException("Error: Main score '%s' is not present in pepXML." % main_score)

	psms = psms.rename(index=str, columns={main_score: 'main_' + main_score})

	# Generate spectrum dataframe
	click.echo("Info: Processing spectra from file %s." % mzxmlfile)
	peaks = read_mzxml(mzxmlfile, psms['scan_id'].unique().tolist())

	return psms, peaks
	