import itertools
import pathlib
from .util import timestamped_echo
import numpy as np
import pandas as pd
import os
import posixpath, ntpath
from statistics import median_low
import click
import re
import time

# Unimod parsing
import xml.etree.cElementTree as ET
from xml.etree.cElementTree import iterparse

# mzXML parsing
import pyopenms as po


class psmtsv:
	relevant_psm_columns = ["Spectrum",
							"Spectrum File",
							"Peptide",
							"Charge",
							"Retention",
							"Delta Mass",
							"Assigned Modifications",
							"Hyperscore",
							"Nextscore",
							"Expectation",
							"Probability",
							"Ion Mobility",
							"Protein",
							"Protein ID",
							"Gene",
							"Mapped Proteins",
							"Mapped Genes",
							"Glycan q-value",
							]
	rank_pattern = re.compile(r"_rank([\d]+)\.pep\.xml")

	def __init__(self, psmtsv_file, unimod, base_name, exclude_range, enable_unannotated, ignore_unannotated, enable_massdiff, decoy_prefix, labile_mods, max_glycan_q):
		self.psmtsv_file = psmtsv_file
		self.base_name = base_name
		self.decoy_prefix = decoy_prefix
		self.labile_mods = labile_mods
		self.max_glycan_q = max_glycan_q
		self.psms = self.parse_psmtsv()
		self.exclude_range = exclude_range
		self.enable_unannotated = enable_unannotated
		self.ignore_unannotated = ignore_unannotated
		self.enable_massdiff = enable_massdiff
		self.match_unimod(unimod)

	def get(self):
		return(self.psms)

	def parse_psmtsv(self):
		# read relevant PSM table columns
		psms = pd.read_csv(self.psmtsv_file, index_col=False, sep='\t', usecols=lambda col: col in set(psmtsv.relevant_psm_columns))

		psms = psms.apply(self.parse_psm_info, axis=1)
		if self.max_glycan_q != 1:
			if "Glycan q-value" in psms.columns:
				psms = psms[(psms['Glycan q-value'] < self.max_glycan_q) | (psms['Glycan q-value'].isnull())]		# include null to include non-glyco peptides
			else:
				timestamped_echo("Warning: glycan q-value filtering requested, but PSM table does not contain the Glycan q-value column. No filtering performed")
		psms = psms.rename(columns={'Charge': 'precursor_charge',
									'Retention': 'retention_time',
									'Delta Mass': 'massdiff',
									'Peptide': 'peptide_sequence',
									'Ion Mobility': 'ion_mobility',
									'Hyperscore': 'var_hyperscore',
									'Nextscore': 'var_nextscore',
									'Expectation': 'var_expect',
									'Probability': 'pep',
									'Spectrum File': 'spectrum_file'
									})
		psms['pep'] = 1 - psms['pep']
		if 'ion_mobility' not in psms:
			psms['ion_mobility'] = np.nan
		if 'hit_rank' not in psms:
			psms['hit_rank'] = 1
		psms = psms.drop(columns=['Spectrum', 'Assigned Modifications', 'Protein', 'Gene', 'Mapped Proteins', 'Mapped Genes', 'Protein ID'])
		if 'Glycan q-value' in psms.columns:
			psms.drop(columns=['Glycan q-value'])
		return psms

	def match_unimod(self, unimod):
		"""
		Match modifications to Unimod as in pepxml. Supports labile modifications
		"""
		def match_modifications(peptide, um, modified_peptide_attribute):
			modified_peptide = peptide['peptide_sequence']

			# parse terminal modifications
			nterm_modification = ""
			if peptide['nterm_modification'] != '':
				nterm_modification = peptide['nterm_modification']
			cterm_modification = ""
			if peptide['cterm_modification'] != '':
				cterm_modification = peptide['cterm_modification']

			# parse closed modifications
			modifications = {}
			if "M|" in peptide[modified_peptide_attribute]:
				for modification in peptide[modified_peptide_attribute].split('|')[1:]:
					site, mass = modification.split('$')
					delta_mass = float(mass)		# psm.tsv mod masses are just the mod, does not include residue mass
					modifications[int(site)] = delta_mass

			massdiff = float(peptide['massdiff'])
			if self.enable_massdiff and (massdiff < self.exclude_range[0] or massdiff > self.exclude_range[1]):
				# parse open modifications
				oms_sequence = peptide['peptide_sequence']
				for site in modifications.keys():
					oms_sequence = oms_sequence[:site-1] + "_" + oms_sequence[site:]

				oms_modifications, nterm_modification, cterm_modification = um.get_oms_id(oms_sequence, peptide['massdiff'], nterm_modification, cterm_modification)
				modifications = {**modifications, **oms_modifications}

			peptide_sequence = peptide['peptide_sequence']
			peptide_length = len(peptide_sequence)
			for site in sorted(modifications, reverse=True):
				positions = ('Anywhere', 'Any N-term', 'Protein N-term') if site == 1 else \
					('Anywhere', 'Any C-term', 'Protein C-term') if site == peptide_length else \
						'Anywhere'
				record_id0 = um.get_id(peptide_sequence[site - 1], positions, modifications[site])
				if isinstance(record_id0, tuple):
					record_id, position = record_id0
				else:
					record_id = record_id0
				is_N_term = isinstance(record_id0, tuple) and position in ('Any N-term', 'Protein N-term')
				if record_id == -1:
					if self.ignore_unannotated:
						return ''		# set empty modified peptide to ignore with later filtering
					elif self.enable_unannotated:
						modified_peptide = "[" + ("+" if modifications[site] > 0 else "") + str(round(modifications[site], 6)) + "]" + modified_peptide \
						if is_N_term else \
						modified_peptide[:site] + "[" + ("+" if modifications[site] > 0 else "") + str(round(modifications[site], 6)) + "]" + modified_peptide[site:]
					else:
						raise click.ClickException("Error: Could not annotate site %s (%s) from peptide %s with delta mass %s." % (site, peptide['peptide_sequence'][site-1], peptide['peptide_sequence'], modifications[site]))
				else:
					modified_peptide = "(UniMod:" + str(record_id) + ")" + modified_peptide \
						if is_N_term else \
						modified_peptide[:site] + "(UniMod:" + str(record_id) + ")" + modified_peptide[site:]

			if nterm_modification != '':
				record_id_nterm = um.get_id("N-term", 'Any N-term', nterm_modification)
				if record_id_nterm == -1:
					record_id_nterm = um.get_id("N-term", 'Protein N-term', nterm_modification)

				if record_id_nterm == -1:
					if self.ignore_unannotated:
						return ''		# set empty modified peptide to ignore with later filtering
					elif self.enable_unannotated:
						modified_peptide = f'.[{round(nterm_modification, 6)}]{modified_peptide}'
					else:
						raise click.ClickException("Error: Could not annotate N-terminus from peptide %s with delta mass %s." % (peptide['peptide_sequence'], nterm_modification))
				else:
					modified_peptide = ".(UniMod:" + str(record_id_nterm) + ")" + modified_peptide

			if cterm_modification != '':
				record_id_cterm = um.get_id("C-term", 'Any C-term', cterm_modification)
				if record_id_cterm == -1:
					record_id_cterm = um.get_id("C-term", 'Protein C-term', cterm_modification)

				if record_id_cterm == -1:
					if self.ignore_unannotated:
						return ''		# set empty modified peptide to ignore with later filtering
					elif self.enable_unannotated:
						modified_peptide = f'{modified_peptide}.[{round(cterm_modification, 6)}]'
					else:
						raise click.ClickException("Error: Could not annotate C-terminus from peptide %s with delta mass %s." % (peptide['peptide_sequence'], cterm_modification))
				else:
					modified_peptide = modified_peptide + ".(UniMod:" + str(record_id_cterm) + ")"

			return modified_peptide

		if self.psms.shape[0] > 0:
			self.psms['modified_peptide'] = self.psms[['peptide_sequence','modifications','nterm_modification','cterm_modification','massdiff']].apply(lambda x: match_modifications(x, unimod, 'modifications'), axis=1)
			self.psms['labile_modified_peptide'] = self.psms[['peptide_sequence','nonlabile_modifications','nterm_modification','cterm_modification','massdiff']].apply(lambda x: match_modifications(x, unimod, 'nonlabile_modifications'), axis=1)
			if self.ignore_unannotated:
				pre_size = len(self.psms)
				self.psms = self.psms[self.psms['modified_peptide'] != '']
				timestamped_echo("Info: Ignored %s PSMs with modifications not matched to Unimod" % (pre_size - len(self.psms)))

	def parse_psm_info(self, psm_series):
		"""
		Perform parsing operations on a PSM entry
		"""
		psm_series = self.parse_spectrum(psm_series)
		psm_series = self.parse_rank(psm_series)
		psm_series = self.parse_assigned_modifications(psm_series)
		psm_series = self.parse_protein_and_gene(psm_series, self.decoy_prefix)
		return psm_series

	def parse_spectrum(self, psm_series):
		splits = psm_series['Spectrum'].split('.')
		psm_series['run_id'] = splits[0]
		psm_series['scan_id'] = int(splits[1])
		return psm_series

	def parse_rank(self, psm_series):
		rank_match = re.search(psmtsv.rank_pattern, psm_series['Spectrum File'])
		if rank_match:
			psm_series['hit_rank'] = int(rank_match.group(1))
		return psm_series

	def parse_assigned_modifications(self, psm_series):
		"""
		Parse the Assigned Modifications column of a psm.tsv table to easyPQP mods.
		Usage: psm_df = psm_df.apply(method, axis=1), where psm_df is the input PSMs dataframe
		Adds the required modification columns to the dataframe
		"""
		modifications = "M"
		nonlabile_modifications = "M"
		nterm_modification = ""
		cterm_modification = ""
		if pd.notnull(psm_series['Assigned Modifications']):
			mods = psm_series['Assigned Modifications'].split(',')
			for mod in mods:
				match = re.search(r"\((-?\d+\.\d+)\)", mod)
				if match:
					mass = match.group(1)
				else:
					click.echo(
						"Error: invalid modification {} in spectrum {} was ignored".format(mod, psm_series['Spectrum']))
					continue

				if 'N-term' in mod:
					nterm_modification = float(mass)
				elif 'C-term' in mod:
					cterm_modification = float(mass)
				else:
					# regular mod
					splits = mod.split('(')
					location = re.search(r"(\d+)", splits[0]).group(1)
					modifications += '|{}${}'.format(location, mass)
					if self.labile_mods != '':
						if self.labile_mods == 'oglyc':
							if not(float(mass) > 140 and psm_series['Peptide'][int(location) - 1] in ['S', 'T']):
								nonlabile_modifications += '|{}${}'.format(location, mass)
						elif self.labile_mods == 'nglyc':
							if not(float(mass) > 140 and psm_series['Peptide'][int(location) - 1] in ['N']):
								nonlabile_modifications += '|{}${}'.format(location, mass)
						elif self.labile_mods == 'nglyc+':
							if not(float(mass) > 140 and psm_series['Peptide'][int(location) - 1] in ['N']):
								nonlabile_modifications += '|{}${}'.format(location, mass)
							else:
								# hard code HexNAc remainder mass as the modification mass for glycan fragment ions rather than full modification mass
								nonlabile_modifications += '|{}${}'.format(location, 203.07937)
					else:
						nonlabile_modifications += '|{}${}'.format(location, mass)
		psm_series['modifications'] = modifications
		psm_series['nonlabile_modifications'] = nonlabile_modifications
		psm_series['nterm_modification'] = nterm_modification
		psm_series['cterm_modification'] = cterm_modification
		return psm_series

	def parse_protein_and_gene(self, psm_series, decoy_prefix):
		"""
		Parse protein and gene IDs from identified and mapped entries and set total_num_proteins and decoy status.
		"""
		psm_series['decoy'] = psm_series['Protein'].startswith(decoy_prefix)

		protein_id = psm_series['Protein ID']
		if pd.notnull(psm_series['Gene']):
			gene_id = psm_series['Gene']
		else:
			gene_id = ''		# contaminants may have empty Gene entry. Ensure string format
		num_total_proteins = 1

		# if mapped proteins/genes not null, capture all entries
		if pd.notnull(psm_series['Mapped Proteins']):
			splits = psm_series['Mapped Proteins'].split(',')
			for split in splits:
				num_total_proteins += 1
				# If UniProt-like IDs, take ID portion only. Otherwise, use whole ID
				uniprot_splits = split.split('|')
				if len(uniprot_splits) > 1:
					protein_id += ';{}'.format(uniprot_splits[1])
				else:
					protein_id += ';{}'.format(split)
		if pd.notnull(psm_series['Mapped Genes']):
			splits = psm_series['Mapped Genes'].split(',')
			for split in splits:
				gene_id += ';{}'.format(split)

		psm_series['protein_id'] = protein_id
		psm_series['gene_id'] = gene_id
		return psm_series


class pepxml:
	def __init__(self, pepxml_file, unimod, base_name, exclude_range, enable_unannotated, enable_massdiff):
		self.pepxml_file = pepxml_file
		self.base_name = base_name
		self.psms = self.parse_pepxml()
		self.exclude_range = exclude_range
		self.enable_unannotated = enable_unannotated
		self.enable_massdiff = enable_massdiff
		self.match_unimod(unimod)

	def get(self):
		return(self.psms)

	def match_unimod(self, unimod):
		def match_modifications(um, peptide):
			monomeric_masses = {"A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694, "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146, "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496, "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203, "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841,
								'U': 150.95363, 'O': 237.14773, 'B': 0, 'J': 0, 'X': 0, 'Z': 0}
			modified_peptide = peptide['peptide_sequence']

			# parse terminal modifications
			nterm_modification = ""
			if peptide['nterm_modification'] != '':
				nterm_modification = peptide['nterm_modification'] - 1.0078
			cterm_modification = ""
			if peptide['cterm_modification'] != '':
				cterm_modification = peptide['cterm_modification'] - 18.0153

			# parse closed modifications
			modifications = {}
			if "M|" in peptide['modifications']:
				for modification in peptide['modifications'].split('|')[1:]:
					site, mass = modification.split('$')
					delta_mass = float(mass) - monomeric_masses[peptide['peptide_sequence'][int(site)-1]]
					modifications[int(site)] = delta_mass

			massdiff = float(peptide['massdiff'])
			if self.enable_massdiff and (massdiff < self.exclude_range[0] or massdiff > self.exclude_range[1]):
				# parse open modifications
				oms_sequence = peptide['peptide_sequence']
				for site in modifications.keys():
					oms_sequence = oms_sequence[:site-1] + "_" + oms_sequence[site:]

				oms_modifications, nterm_modification, cterm_modification = um.get_oms_id(oms_sequence, peptide['massdiff'], nterm_modification, cterm_modification)
				modifications = {**modifications, **oms_modifications}

			peptide_sequence = peptide['peptide_sequence']
			peptide_length = len(peptide_sequence)
			for site in sorted(modifications, reverse=True):
				positions = ('Anywhere', 'Any N-term', 'Protein N-term') if site == 1 else \
					('Anywhere', 'Any C-term', 'Protein C-term') if site == peptide_length else \
						'Anywhere'
				record_id0 = um.get_id(peptide_sequence[site - 1], positions, modifications[site])
				if isinstance(record_id0, tuple):
					record_id, position = record_id0
				else:
					record_id = record_id0
				is_N_term = isinstance(record_id0, tuple) and position in ('Any N-term', 'Protein N-term')
				if record_id == -1:
					if self.enable_unannotated:
						modified_peptide = "[" + ("+" if modifications[site] > 0 else "") + str(round(modifications[site], 6)) + "]" + modified_peptide \
						if is_N_term else \
						modified_peptide[:site] + "[" + ("+" if modifications[site] > 0 else "") + str(round(modifications[site], 6)) + "]" + modified_peptide[site:]
					else:
						raise click.ClickException("Error: Could not annotate site %s (%s) from peptide %s with delta mass %s." % (site, peptide['peptide_sequence'][site-1], peptide['peptide_sequence'], modifications[site]))
				else:
					modified_peptide = "(UniMod:" + str(record_id) + ")" + modified_peptide \
						if is_N_term else \
						modified_peptide[:site] + "(UniMod:" + str(record_id) + ")" + modified_peptide[site:]

			if nterm_modification != '':
				record_id_nterm = um.get_id("N-term", 'Any N-term', nterm_modification)
				if record_id_nterm == -1:
					record_id_nterm = um.get_id("N-term", 'Protein N-term', nterm_modification)

				if record_id_nterm == -1:
					if self.enable_unannotated:
						modified_peptide = f'.[{round(nterm_modification, 6)}]{modified_peptide}'
					else:
						raise click.ClickException("Error: Could not annotate N-terminus from peptide %s with delta mass %s." % (peptide['peptide_sequence'], nterm_modification))
				else:
					modified_peptide = ".(UniMod:" + str(record_id_nterm) + ")" + modified_peptide

			if cterm_modification != '':
				record_id_cterm = um.get_id("C-term", 'Any C-term', cterm_modification)
				if record_id_cterm == -1:
					record_id_cterm = um.get_id("C-term", 'Protein C-term', cterm_modification)

				if record_id_cterm == -1:
					if self.enable_unannotated:
						modified_peptide = f'{modified_peptide}.[{round(cterm_modification, 6)}]'
					else:
						raise click.ClickException("Error: Could not annotate C-terminus from peptide %s with delta mass %s." % (peptide['peptide_sequence'], cterm_modification))
				else:
					modified_peptide = modified_peptide + ".(UniMod:" + str(record_id_cterm) + ")"


			return modified_peptide

		if self.psms.shape[0] > 0:
			self.psms['modified_peptide'] = self.psms[['peptide_sequence','modifications','nterm_modification','cterm_modification','massdiff']].apply(lambda x: match_modifications(unimod, x), axis=1)

	def parse_pepxml(self):
		peptides = []
		namespaces = {'pepxml_ns': "http://regis-web.systemsbiology.net/pepXML"}
		ET.register_namespace('', "http://regis-web.systemsbiology.net/pepXML")

		context = iterparse(self.pepxml_file, events=("end",))

		for event, elem in context:
			if elem.tag == "{http://regis-web.systemsbiology.net/pepXML}msms_run_summary":
				base_name = basename_spectralfile(posixpath.basename(ntpath.basename(elem.attrib['base_name'])))

				# only proceed if base_name matches
				if base_name == self.base_name:
					# find decoy prefix
					decoy_prefix = ""
					for search_summary in elem.findall('.//pepxml_ns:search_summary', namespaces):
						for parameter in search_summary.findall('.//pepxml_ns:parameter', namespaces):
							if parameter.attrib['name'] == 'decoy_prefix':
								decoy_prefix = parameter.attrib['value']

					# go through all spectrum queries
					for spectrum_query in elem.findall('.//pepxml_ns:spectrum_query', namespaces):
						index = spectrum_query.attrib['index']
						start_scan = spectrum_query.attrib['start_scan']
						end_scan = spectrum_query.attrib['end_scan']
						assumed_charge = spectrum_query.attrib['assumed_charge']
						retention_time_sec = spectrum_query.attrib['retention_time_sec']

						ion_mobility = np.nan
						if 'ion_mobility' in spectrum_query.attrib:
							ion_mobility = spectrum_query.attrib['ion_mobility']

						for search_result in spectrum_query.findall(".//pepxml_ns:search_result", namespaces):
							prev_pep = None
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
										proteins[prot.split("|")[1]] = prot.split("|")[2].split(" ")[0].split("_")[0]
									else:
										proteins[prot] = prot

								protein = ""
								gene = ""

								for key in sorted(proteins):
									if protein == "":
										protein = key
									else:
										protein = protein + ";" + key

									if gene == "":
										gene = proteins[key]
									else:
										gene = gene + ";" + proteins[key]

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
											prev_pep = scores["pep"]
									elif analysis_result.attrib['analysis'] == 'peptideprophet':
										for peptideprophet_result in analysis_result.findall('.//pepxml_ns:peptideprophet_result', namespaces):
											scores["pep"] = 1.0 - float(peptideprophet_result.attrib['probability'])
											prev_pep = scores["pep"]

								if prev_pep is not None and "pep" not in scores:
									# If 2 search hits have the same rank only the first one has the analysis_result explicitly written out.
									scores["pep"] = prev_pep

								if prev_pep is not None:
									peptides.append({**{'run_id': base_name, 'scan_id': int(start_scan), 'hit_rank': int(hit_rank), 'massdiff': float(massdiff), 'precursor_charge': int(assumed_charge), 'retention_time': float(retention_time_sec), 'ion_mobility': float(ion_mobility), 'peptide_sequence': peptide, 'modifications': modifications, 'nterm_modification': nterm_modification, 'cterm_modification': cterm_modification, 'protein_id': protein, 'gene_id': gene, 'num_tot_proteins': num_tot_proteins, 'decoy': is_decoy}, **scores})
				elem.clear()

		df = pd.DataFrame(peptides)

		return(df)

class idxml:
    def __init__(self, idxml_file, base_name):
        self.idxml_file = idxml_file
        self.base_name = base_name
        self.psms = self.parse_idxml()
        #self.exclude_range = exclude_range
        #self.match_unimod(unimod)

    def get(self):
        return(self.psms)

    def parse_idxml(self):
        peptides = []
        proteins = []
        scores = {}

        parsed_peptides = []

        po.IdXMLFile().load(self.idxml_file, proteins, peptides)

        for i, p in enumerate(peptides):

            #percolator probability
            scores["q_value"] = float(p.getHits()[0].getScore() if p.getScoreType() == "q-value" else 0)
            scores["pep"] = float(p.getHits()[0].getScore() if p.getScoreType() == "PEP" else 0)

            try:
                parsed_peptides.append({**{'run_id': self.base_name,
                                            'scan_id': int(get_scan(p.getMetaValue("spectrum_reference"), i)),
                                            'hit_rank': int(p.getHits()[0].getRank()),
                                            'massdiff': float(0),
                                            'precursor_charge': int(p.getHits()[0].getCharge()),
                                            'retention_time': float(p.getRT()),
                                            'ion_mobility': float(p.getMetaValue("IM")) if p.getMetaValue("IM") is not None else np.nan,
                                            'modified_peptide': p.getHits()[0].getSequence().toUniModString().decode("utf-8"),
                                            'peptide_sequence': p.getHits()[0].getSequence().toUnmodifiedString().decode("utf-8"),
                                            'modifications': '-',
                                            'nterm_modification': '-',
                                            'cterm_modification': '-',
                                            'protein_id': ','.join([prot.getProteinAccession().decode("utf-8") for prot in p.getHits()[0].getPeptideEvidences()]),
                                            'gene_id': '-',
                                            'num_tot_proteins': len([prot.getProteinAccession() for prot in p.getHits()[0].getPeptideEvidences()]),
                                            'decoy': p.getHits()[0].getMetaValue('target_decoy').decode("utf-8")=='decoy'}, **scores})

            except AttributeError:
                parsed_peptides.append({**{'run_id': self.base_name,
                                            'scan_id': int(get_scan(p.getMetaValue("spectrum_reference"), i)),
                                            'hit_rank': int(p.getHits()[0].getRank()),
                                            'massdiff': float(0),
                                            'precursor_charge': int(p.getHits()[0].getCharge()),
                                            'retention_time': float(p.getRT()),
                                            'ion_mobility': float(p.getMetaValue("IM")) if p.getMetaValue("IM") is not None else np.nan,
                                            'modified_peptide': p.getHits()[0].getSequence().toUniModString(),
                                            'peptide_sequence': p.getHits()[0].getSequence().toUnmodifiedString(),
                                            'modifications': '-',
                                            'nterm_modification': '-',
                                            'cterm_modification': '-',
                                            'protein_id': ','.join([prot.getProteinAccession() for prot in p.getHits()[0].getPeptideEvidences()]),
                                            'gene_id': '-',
                                            'num_tot_proteins': len([prot.getProteinAccession() for prot in p.getHits()[0].getPeptideEvidences()]),
                                            'decoy': p.getHits()[0].getMetaValue('target_decoy')=='decoy'}, **scores})

        df = pd.DataFrame(parsed_peptides)

        return (df)

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
		sites = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V','N-term','C-term', 'B', 'J', 'X', 'Z']
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
		ptms_site = self.ptms[site]
		search_multiple_positions = isinstance(position, (list, tuple))

		kvs = itertools.chain.from_iterable((((k, p), v) for k,v in ptms_site[p].items()) for p in position) \
			if search_multiple_positions else \
			ptms_site[position].items()
		for key, value in kvs:
			delta_mod = abs(value - float(delta_mass))
			if delta_mod < self.max_delta:
				candidates_key = candidates.setdefault(key, delta_mod)
				if delta_mod < candidates_key:
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


good_patterns = [re.compile('controllerType=0 controllerNumber=1 scan=([0-9]+)'),
				re.compile('frame=([0-9]+)'),
				re.compile('function=[0-9]+ process=[0-9]+ scan=([0-9]+)'),
				re.compile('jobRun=[0-9]+ spotLabel=[^ ]+ spectrum=([0-9]+)'),
				re.compile('([0-9]+)'),
				re.compile('scan=([0-9]+)'),
				re.compile('spectrum=([0-9]+)'),
				re.compile('scanId=([0-9]+)'),
				re.compile('index=([0-9]+)')]

bad_patterns = [re.compile('controllerType=[1-9] controllerNumber=1 scan=[0-9]+'),
				re.compile('controllerType=[0-9] controllerNumber=1 scan=[0-9]+ demux=[0-9]+')]


def get_scan(e: str, fallback_num: int):
	if e is None or e == "":
		return fallback_num

	for bad_pattern in bad_patterns:
		if bad_pattern.fullmatch(e) is not None:
			return fallback_num

	for good_pattern in good_patterns:
		x = good_pattern.fullmatch(e)
		if x is not None:
			return int(x.group(1))

	return fallback_num


def read_mzml_or_mzxml_impl(input_map, psms, theoretical, max_delta_ppm, filetype):
	assert filetype in ('mzml', 'mzxml')

	timestamped_echo("Info: Collecting PSMs...")

	input_map = {get_scan(getNativeID, idx + 1): e for idx, (getNativeID, e) in enumerate(input_map)}
	import concurrent.futures
	nthreads = min(os.cpu_count(), 5)
	def f(psms):
		peaks_list = []
		for scan_id, modified_peptide, precursor_charge in psms.itertuples(index=None):
			peaks_list.append(
				psm_df(input_map, theoretical, max_delta_ppm, scan_id, modified_peptide, precursor_charge))
		return peaks_list
	with concurrent.futures.ThreadPoolExecutor(nthreads) as exe:
		l = exe.map(f, np.array_split(psms, nthreads))
	peaks_list = sum(l, [])

	timestamped_echo("Info: Got %d PSMs" % len(peaks_list))
	timestamped_echo("Info: Generating Transitions...")

	if len(peaks_list) > 0:
		reps = np.array([e[0] for e in peaks_list])
		transitions = pd.DataFrame({'fragment': np.concatenate([e[1] for e in peaks_list]),
									'product_mz': np.concatenate([e[2] for e in peaks_list]),
									'intensity': np.concatenate([e[3] for e in peaks_list]),
									'scan_id': np.repeat([e[4] for e in peaks_list], reps),
									'precursor_mz': np.repeat([e[5] for e in peaks_list], reps),
									'modified_peptide': np.repeat([e[6] for e in peaks_list], reps),
									'precursor_charge': np.repeat([e[7] for e in peaks_list], reps)})
		# Multiple peaks might be identically annotated, only use most intense
		transitions = transitions.groupby(['scan_id','modified_peptide','precursor_charge','precursor_mz','fragment','product_mz'])['intensity'].max().reset_index()
	else:
		transitions = pd.DataFrame({'scan_id': [], 'modified_peptide': [], 'precursor_charge': [], 'precursor_mz': [], 'fragment': [], 'product_mz': [], 'intensity': []})

	timestamped_echo("Info: Generated %d transitions" % len(transitions))

	return(transitions)


def read_mgf(tims_mgf_path, psms, theoretical, max_delta_ppm):
	'''
	read TimsTOF's MGF or MSFragger's calibrated.mgf
	'''
	is_tims = not tims_mgf_path.endswith('calibrated.mgf')
	import mmap
	record_pattern = re.compile(b'''BEGIN IONS\r?
(.*?)
END IONS''', re.MULTILINE | re.DOTALL)
	scan_num_pattern_Bruker = re.compile(rb'TITLE=Cmpd\s+([0-9]+),')
	scan_num_pattern_general = re.compile(rb'TITLE=.+?\.([0-9]+?)\.')
	peaks_pattern = re.compile(rb'^([\d.]+)\s+([\d.]+)', re.MULTILINE)

	tims_data = {}
	with open(tims_mgf_path, "rb") as f:
		mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
		for e in record_pattern.finditer(mm):
			rec = e.group(1)
			scan_num_findall = (scan_num_pattern_Bruker.findall(rec) + scan_num_pattern_general.findall(rec))
			if len(scan_num_findall) in [1, 2]:
				scan_num = int(scan_num_findall[0])
			else:
				raise RuntimeError("Cannot find Cmpd number from " + rec.decode()
								   if is_tims else
								   "Cannot find scan number from " + rec.decode())
			tims_data[scan_num] = np.array(peaks_pattern.findall(rec), dtype=float)

	peaks_list = []
	for scan_id, modified_peptide, precursor_charge in psms.itertuples(index=False):
		peaks_list.append(psm_df_mgf(tims_data, theoretical, max_delta_ppm, scan_id, modified_peptide, precursor_charge))

	if len(peaks_list) > 0:
		reps = np.array([e[0] for e in peaks_list])
		transitions = pd.DataFrame({'fragment': np.concatenate([e[1] for e in peaks_list]),
									'product_mz': np.concatenate([e[2] for e in peaks_list]),
									'intensity': np.concatenate([e[3] for e in peaks_list]),
									'scan_id': np.repeat([e[4] for e in peaks_list], reps),
									'precursor_mz': np.repeat([e[5] for e in peaks_list], reps),
									'modified_peptide': np.repeat([e[6] for e in peaks_list], reps),
									'precursor_charge': np.repeat([e[7] for e in peaks_list], reps)})
		# Multiple peaks might be identically annotated, only use most intense
		transitions = transitions.groupby(['scan_id','modified_peptide','precursor_charge','precursor_mz','fragment','product_mz'])['intensity'].max().reset_index()
	else:
		transitions = pd.DataFrame({'scan_id': [], 'modified_peptide': [], 'precursor_charge': [], 'precursor_mz': [], 'fragment': [], 'product_mz': [], 'intensity': []})
	return transitions


def annotate_mass(mass, ionseries, max_delta_ppm):
	top_delta = 30
	ions, ion_masses = ionseries
	ppms = np.abs((mass - ion_masses) / ion_masses * 1e6)
	idx = ppms.argmin()
	if ppms[idx] < min(max_delta_ppm, top_delta):
		return ions[idx], ion_masses[idx]
	return None, None


def psm_df(input_map, theoretical, max_delta_ppm, scan_id, modified_peptide, precursor_charge):
	ionseries = theoretical[modified_peptide][precursor_charge]

	spectrum = input_map[scan_id]

	fragments, product_mzs, intensities = annotate_mass_spectrum_numba(ionseries, max_delta_ppm, spectrum)
	# Baseline normalization to highest annotated peak
	max_intensity = np.amax(intensities, initial=0.0)
	if max_intensity > 0:
		intensities /= max_intensity
		intensities *= 10000
	return [len(fragments), fragments, product_mzs, intensities,
			scan_id,
			po.AASequence.fromString(po.String(modified_peptide))
				.getMonoWeight(po.Residue.ResidueType.Full, precursor_charge) / precursor_charge,
			modified_peptide, precursor_charge]

def psm_df_mgf(input_map, theoretical, max_delta_ppm, scan_id, modified_peptide, precursor_charge):
	ionseries = theoretical[modified_peptide][precursor_charge]

	spectrum = input_map[scan_id]

	top_delta = 30
	ions, ion_masses = ionseries

	mzs0, intensities0 = spectrum[:, 0], spectrum[:, 1]
	ppms = np.abs((mzs0[:, np.newaxis] - ion_masses) / ion_masses * 1e6)
	idx_mask = (ppms < min(max_delta_ppm, top_delta)).any(1)
	idx = ppms[idx_mask].argmin(1)
	fragments, product_mzs, intensities = ions[idx], ion_masses[idx], intensities0[idx_mask]

	# Baseline normalization to highest annotated peak
	max_intensity = np.amax(intensities, initial=0.0)
	if max_intensity > 0:
		intensities /= max_intensity
		intensities *= 10000
	return [len(fragments), fragments, product_mzs, intensities,
			scan_id,
			po.AASequence.fromString(po.String(modified_peptide))
				.getMonoWeight(po.Residue.ResidueType.Full, precursor_charge) / precursor_charge,
			modified_peptide, precursor_charge]


def annotate_mass_spectrum(ionseries, max_delta_ppm, spectrum):
	top_delta = 30
	ions, ion_masses = ionseries

	mzs0, intensities0 = spectrum#.get_peaks()
	ppms = np.abs((mzs0[:, np.newaxis] - ion_masses) / ion_masses * 1e6)
	idx_mask = (ppms < min(max_delta_ppm, top_delta)).any(1)
	idx = ppms[idx_mask].argmin(1)
	return ions[idx], ion_masses[idx], intensities0[idx_mask]

import numba
@numba.jit(nopython=True, nogil=True, fastmath=True)
def annotate_mass_spectrum_numba(ionseries, max_delta_ppm, spectrum):
	top_delta = min(max_delta_ppm, 30)
	ions, ion_masses = ionseries
	mzs0, intensities0 = spectrum
	idx_ions = np.empty_like(intensities0, dtype=np.uint32)
	idx_peaks = np.zeros_like(intensities0, dtype=np.bool_)
	for si, mz in enumerate(mzs0):
		min_ion_idx = -1
		min_ppm = top_delta
		for ii, ion_mass in enumerate(ion_masses):
			ppm = np.abs(mz - ion_mass) / ion_mass * 1e6
			if ppm < min_ppm:
				min_ppm = ppm
				min_ion_idx = ii
		if min_ion_idx != -1:
			idx_ions[si] = min_ion_idx
			idx_peaks[si] = True
	return ions[idx_ions[idx_peaks]], ion_masses[idx_ions[idx_peaks]], intensities0[idx_peaks]


def generate_ionseries(peptide_sequence, precursor_charge, fragment_charges=[1,2,3,4], fragment_types=['b','y'], enable_specific_losses = False, enable_unspecific_losses = False, precision_digits = 6):
	peptide = po.AASequence.fromString(po.String(peptide_sequence))
	sequence = peptide.toUnmodifiedString()

	# this dict is compatible with DIA-NN
	unspecific_losses = dict()
	unspecific_losses["H2O"] = 18.0106
	unspecific_losses["NH3"] = 17.0265
	# unspecific_losses["CO"] = 27.9949
	# unspecific_losses["H3PO4"] = 97.9769
	# unspecific_losses["H4COS"] = 63.9983

	fragments = {}

	for fragment_type in fragment_types:
		for fragment_charge in fragment_charges:
			if fragment_charge <= precursor_charge:
				for fragment_ordinal in range(1,len(sequence)):
					mass = ion = None
					if fragment_type == 'a':
						ion = peptide.getPrefix(fragment_ordinal)
						mass = ion.getMonoWeight(po.Residue.ResidueType.AIon, fragment_charge) / fragment_charge;
					elif fragment_type == 'b':
						ion = peptide.getPrefix(fragment_ordinal)
						mass = ion.getMonoWeight(po.Residue.ResidueType.BIon, fragment_charge) / fragment_charge;
					elif fragment_type == 'c':
						ion = peptide.getPrefix(fragment_ordinal)
						mass = ion.getMonoWeight(po.Residue.ResidueType.CIon, fragment_charge) / fragment_charge;
					elif fragment_type == 'x':
						ion = peptide.getSuffix(fragment_ordinal)
						mass = ion.getMonoWeight(po.Residue.ResidueType.XIon, fragment_charge) / fragment_charge;
					elif fragment_type == 'y':
						ion = peptide.getSuffix(fragment_ordinal)
						mass = ion.getMonoWeight(po.Residue.ResidueType.YIon, fragment_charge) / fragment_charge;
					elif fragment_type == 'z':
						ion = peptide.getSuffix(fragment_ordinal)
						mass = ion.getMonoWeight(po.Residue.ResidueType.ZIon, fragment_charge) / fragment_charge;
					else:
						raise RuntimeError(f'fragment type "{fragment_type}" is not in (a,b,c,x,y,z)')
					
					# Workaround
					# If two fragment ions have identical product m/z, this can lead to annotation inconsistencies.
					# E.g. .(UniMod:1)ADQLTEEQIAEFK+2 and the corresponding b5^1 and b10^2 ions.
					# We thus first generate a dict with product m/z as index and the reserve the dict.
					# Then do not update the dicts if the product m/z is already present.
					# This leads to only the lower charged product to be annotated and reported.

					# Standard fragment ions
					k = round(mass, precision_digits)
					if k not in fragments:
						fragments[k] = fragment_type + str(fragment_ordinal) + "^" + str(fragment_charge)

					# unspecific losses that are compatible with DIA-NN
					if enable_unspecific_losses:
						for loss in unspecific_losses:
							k = round(mass - (unspecific_losses[loss] / fragment_charge), precision_digits)
							if k not in fragments:
								fragments[k] = fragment_type + str(fragment_ordinal) + "-" + loss + "^" + str(fragment_charge)

					# specific losses that are hardcoded in OpenMS
					if enable_specific_losses:
						for lossfragment_ordinal in range(1, ion.size()):
							if ion.getResidue(lossfragment_ordinal).hasNeutralLoss():
								losses = ion.getResidue(lossfragment_ordinal).getLossFormulas()
								for loss in losses:
									loss_type = loss.toString()
									if loss_type not in unspecific_losses:
										k = round(mass - (loss.getMonoWeight() / fragment_charge), precision_digits)
										fragments[k] = fragment_type + str(fragment_ordinal) + "-" + loss_type + "^" + str(fragment_charge)

	# flip key - values for workaround
	fragments = {value: key for key, value in fragments.items()}
	return np.array(list(fragments.keys())), np.fromiter(fragments.values(), float, len(fragments))

def parse_pepxmls(pepxmlfile_list, um, base_name, exclude_range, enable_unannotated, enable_massdiff, fragment_charges, fragment_types, enable_specific_losses, enable_unspecific_losses, precision_digits):
	psmslist = []
	for pepxmlfile in pepxmlfile_list:
		if pepxmlfile.casefold().endswith(('.pepxml', '.pep.xml')):
			timestamped_echo(f"Info: Parsing pepXML: {pepxmlfile}")
			px = pepxml(pepxmlfile, um, base_name, exclude_range, enable_unannotated, enable_massdiff)
		elif pepxmlfile.lower().endswith('idxml'):
			timestamped_echo(f"Info: Parsing idXML: {pepxmlfile}")
			px = idxml(pepxmlfile, base_name)
		else:
			timestamped_echo('unknown format of pepxml identification file')
		psms = px.get()
		if psms.shape[0] > 0:
			rank = re.compile(r'_rank([0-9]+)\.').search(pathlib.Path(pepxmlfile).name)
			rank_str = '' if rank is None else '_rank' + rank.group(1)
			psms['group_id'] = psms['run_id'] + "_" + psms['scan_id'].astype(str) + rank_str
			psmslist.append(psms)
		timestamped_echo(f"Info: {psms.shape[0]} PSMs parsed.")

	out_df = pd.DataFrame()
	theoretical = None
	if len(psmslist) > 0:
		out_df = pd.concat(psmslist)
		if out_df.shape[0] > 0:
			timestamped_echo("Info: Generate theoretical spectra.")
			theoretical = {}
			for modified_peptide, precursor_charge in out_df[['modified_peptide','precursor_charge']].drop_duplicates().itertuples(index=False):
				theoretical.setdefault(modified_peptide, {})[precursor_charge] = generate_ionseries(modified_peptide, precursor_charge, fragment_charges, fragment_types, enable_specific_losses, enable_unspecific_losses, precision_digits)
	return out_df, theoretical

class MSCallback:
	def __init__(self):
		self.id_peaks_map = []

	def setExperimentalSettings(self, s):
		pass

	def setExpectedSize(self, a, b):
		pass

	def consumeChromatogram(self, c):
		pass
	def consumeSpectrum(self, s):
		self.id_peaks_map.append((s.getNativeID(), s.get_peaks()))

def parse_psms(psm_file_list, um, base_name, exclude_range, enable_unannotated, ignore_unannotated, enable_massdiff, fragment_charges, fragment_types, enable_specific_losses, enable_unspecific_losses, decoy_prefix, precision_digits, labile_mods, max_glycan_q):
	"""
	Parsing method with psm.tsv as the primary input instead of pepxml/idxml
	"""
	psmslist = []
	for psmfile in psm_file_list:
		px = psmtsv(psmfile, um, base_name, exclude_range, enable_unannotated, ignore_unannotated, enable_massdiff, decoy_prefix, labile_mods, max_glycan_q)
		psms = px.get()
		if psms['spectrum_file'].str.extract(r'_rank(\d+)').isna().all().all():
			psms['group_id'] = psms['run_id'] + "_" + psms['scan_id'].astype(str)
		else:
			psms['group_id'] = psms['run_id'] + "_" + psms['scan_id'].astype(str) + "_rank" + psms['hit_rank'].astype(str)
		click.echo(f"Info: Done parsing psm.tsv: {psmfile}")
		psmslist.append(psms)
	psms = pd.concat(psmslist)
	theoretical = None
	if psms.shape[0] > 0:
		# Generate theoretical spectra
		click.echo("Info: Generate theoretical spectra.")
		theoretical = {}
		if labile_mods != '':
			for modified_peptide, precursor_charge, labile_peptide in psms[['modified_peptide', 'precursor_charge', 'labile_modified_peptide']].drop_duplicates().itertuples(index=False):
				theoretical.setdefault(modified_peptide, {})[precursor_charge] = generate_ionseries(labile_peptide, precursor_charge, fragment_charges, fragment_types, enable_specific_losses, enable_unspecific_losses, precision_digits)
		else:
			for modified_peptide, precursor_charge in psms[['modified_peptide', 'precursor_charge']].drop_duplicates().itertuples(index=False):
				theoretical.setdefault(modified_peptide, {})[precursor_charge] = generate_ionseries(modified_peptide, precursor_charge, fragment_charges, fragment_types, enable_specific_losses, enable_unspecific_losses, precision_digits)
	return psms, theoretical


def get_map_mzml_or_mzxml(path: str, filetype):
	assert filetype in ('mzml', 'mzxml')
	fh = po.MzMLFile() if filetype=='mzml' else po.MzXMLFile()
	consumer = MSCallback()
	fh.transform(path, consumer)
	return consumer.id_peaks_map

def conversion(pepxmlfile_list, spectralfile, unimodfile, exclude_range, max_delta_unimod, max_delta_ppm, enable_unannotated, enable_massdiff, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, max_psm_pep, precision_digits):
	# Parse basename
	base_name = basename_spectralfile(spectralfile)
	timestamped_echo("Info: Parsing run %s." % base_name)

	# Initialize UniMod
	um = unimod(unimodfile, max_delta_unimod)

	psms, theoretical = parse_pepxmls(pepxmlfile_list, um, base_name, exclude_range, enable_unannotated, enable_massdiff, fragment_charges, fragment_types, enable_specific_losses, enable_unspecific_losses, precision_digits)

	timestamped_echo("Info: Processing spectra from file %s." % spectralfile)

	if spectralfile.lower().endswith(".mzxml"):
		input_map = get_map_mzml_or_mzxml(spectralfile, 'mzxml')
	elif spectralfile.casefold().endswith(".mzml"):
		input_map = get_map_mzml_or_mzxml(spectralfile, 'mzml')
	else:
		input_map = None

	timestamped_echo("Info: Loaded %d spectra" % len(input_map))

	if psms.shape[0] > 0:
		# Generate spectrum dataframe
		psms = psms[psms['pep'] <= max_psm_pep]
		if spectralfile.lower().endswith(".mzxml"):
			peaks = read_mzml_or_mzxml_impl(input_map, psms[['scan_id','modified_peptide','precursor_charge']], theoretical, max_delta_ppm, 'mzxml')
		elif spectralfile.casefold().endswith(".mzml"):
			peaks = read_mzml_or_mzxml_impl(input_map, psms[['scan_id', 'modified_peptide', 'precursor_charge']], theoretical, max_delta_ppm, 'mzml')
		elif spectralfile.lower().endswith(".mgf"):
			peaks = read_mgf(spectralfile, psms[['scan_id', 'modified_peptide', 'precursor_charge']], theoretical, max_delta_ppm)

		# Round floating numbers
		peaks = peaks.round(6)

		# Filter PSMs not having any peaks to avoid the "library" command picking the PSM (with the smallest PEP) without any peaks.
		psms = psms[psms['scan_id'].isin(peaks['scan_id'].unique())]

		return psms, peaks
	else:
		return (pd.DataFrame({'run_id': [],
							  'scan_id': [],
							  'hit_rank': [],
							  'massdiff': [],
							  'precursor_charge': [],
							  'retention_time': [],
							  'ion_mobility': [],
							  'peptide_sequence': [],
							  'protein_id': [],
							  'gene_id': [],
							  'num_tot_proteins': [],
							  'decoy': [],
							  'pep': [],
							  'modified_peptide': [],
							  'group_id': []}),
				pd.DataFrame({'scan_id': [],
							  'modified_peptide': [],
							  'precursor_charge': [],
							  'precursor_mz': [],
							  'fragment': [],
							  'product_mz': [],
							  'intensity': []}))


def conversion_psm(psm_file_list, spectralfile, unimodfile, exclude_range, max_delta_unimod, max_delta_ppm, enable_unannotated, ignore_unannotated, enable_massdiff, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, max_psm_pep, decoy_prefix, precision_digits, labile_mods, max_glycan_q):
	# Parse basename
	base_name = basename_spectralfile(spectralfile)
	timestamped_echo("Info: Parsing run %s." % base_name)

	# Initialize UniMod
	um = unimod(unimodfile, max_delta_unimod)
	import concurrent.futures

	timestamped_echo("Info: Processing spectra from file %s." % spectralfile)

	exe = concurrent.futures.ProcessPoolExecutor(1)
	psms_fut = exe.submit(parse_psms, psm_file_list, um, base_name, exclude_range, enable_unannotated, ignore_unannotated, enable_massdiff, fragment_charges, fragment_types, enable_specific_losses, enable_unspecific_losses, decoy_prefix, precision_digits, labile_mods, max_glycan_q)
	time.sleep(1)  # allow the process to execute first before using pyOpenMS to read files
	if spectralfile.lower().endswith(".mzxml"):
		input_map = get_map_mzml_or_mzxml(spectralfile, 'mzxml')
	elif spectralfile.casefold().endswith(".mzml"):
		input_map = get_map_mzml_or_mzxml(spectralfile, 'mzml')
	else:
		input_map = None

	timestamped_echo("Info: Loaded %d spectra" % len(input_map))
	timestamped_echo("Info: Processing PSM file...")

	# Continue if any PSMS are present
	psms, theoretical = psms_fut.result()
	exe.shutdown()

	if psms.shape[0] > 0:
		# Generate spectrum dataframe
		click.echo("Info: Processing spectra from file %s." % spectralfile)
		psms = psms[psms['pep'] <= max_psm_pep]
		if spectralfile.lower().endswith(".mzxml"):
			peaks = read_mzml_or_mzxml_impl(input_map, psms[['scan_id','modified_peptide','precursor_charge']], theoretical, max_delta_ppm, 'mzxml')
		elif spectralfile.casefold().endswith(".mzml"):
			peaks = read_mzml_or_mzxml_impl(input_map, psms[['scan_id', 'modified_peptide', 'precursor_charge']], theoretical, max_delta_ppm, 'mzml')
		elif spectralfile.lower().endswith(".mgf"):
			peaks = read_mgf(spectralfile, psms[['scan_id', 'modified_peptide', 'precursor_charge']], theoretical, max_delta_ppm)

		# Round floating numbers
		peaks = peaks.round(6)

		# Filter PSMs not having any peaks to avoid the "library" command picking the PSM (with the smallest PEP) without any peaks.
		psms = psms[psms['scan_id'].isin(peaks['scan_id'].unique())]

		return psms, peaks
	else:
		return (pd.DataFrame({'run_id': [],
							 'scan_id': [],
							 'hit_rank': [],
							 'massdiff': [],
							 'precursor_charge': [],
							 'retention_time': [],
							 'ion_mobility': [],
							 'peptide_sequence': [],
							 'protein_id': [],
							 'gene_id': [],
							 'num_tot_proteins': [],
							 'decoy': [],
							 'pep': [],
							 'modified_peptide': [],
							 'group_id': []}),
				pd.DataFrame({'scan_id': [],
							  'modified_peptide': [],
							  'precursor_charge': [],
							  'precursor_mz': [],
							  'fragment': [],
							  'product_mz': [],
							  'intensity': []}))


def drop_psm_columns(psms):
	t = ['run_id',
		 'scan_id',
		 'hit_rank',
		 'massdiff',
		 'precursor_charge',
		 'retention_time',
		 'ion_mobility',
		 'peptide_sequence',
		 'protein_id',
		 'gene_id',
		 'num_tot_proteins',
		 'decoy',
		 'pep',
		 'modified_peptide',
		 'group_id']
	tt = list(set(t) & set(psms.columns))
	return psms.loc[:, tt]


def basename_spectralfile(spectralfile):
	'''
	take the basename of a spectral filename
	and strip trailing `_calibrated` if any, for MSFragger's MGF.
	:param spectralfile: name of spectral file
	:return: basename without trailing `_calibrated`
	'''
	x = os.path.splitext(os.path.basename(spectralfile))[0]
	# get basename without _(un)calibrated suffix, if any
	return re.compile('(.+?)(?:_(?:un)?calibrated)?').fullmatch(x)[1]
