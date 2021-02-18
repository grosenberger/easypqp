import itertools
import pathlib

import numpy as np
import pandas as pd
import os
import posixpath, ntpath
from statistics import median_low
import click
import re

# Unimod parsing
import xml.etree.cElementTree as ET
from xml.etree.cElementTree import iterparse

# mzXML parsing
import pyopenms as po

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
			monomeric_masses = {"A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694, "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146, "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496, "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203, "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841}
			modified_peptide = peptide['peptide_sequence']

			# parse terminal modifications
			nterm_modification = ""
			if peptide['nterm_modification'] is not "":
				nterm_modification = peptide['nterm_modification'] - 1.0078
			cterm_modification = ""
			if peptide['cterm_modification'] is not "":
				cterm_modification = peptide['cterm_modification']

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
						modified_peptide = "[" + str(modifications[site]) + "]" + modified_peptide \
						if is_N_term else \
						modified_peptide[:site] + "[" + str(modifications[site]) + "]" + modified_peptide[site:]
					else:
						raise click.ClickException("Error: Could not annotate site %s (%s) from peptide %s with delta mass %s." % (site, peptide['peptide_sequence'][site-1], peptide['peptide_sequence'], modifications[site]))
				else:
					modified_peptide = "(UniMod:" + str(record_id) + ")" + modified_peptide \
						if is_N_term else \
						modified_peptide[:site] + "(UniMod:" + str(record_id) + ")" + modified_peptide[site:]

			if nterm_modification is not "":
				record_id_nterm = um.get_id("N-term", 'Any N-term', nterm_modification)
				if record_id_nterm == -1:
					record_id_nterm = um.get_id("N-term", 'Protein N-term', nterm_modification)

				if record_id_nterm == -1:
					if self.enable_unannotated:
						modified_peptide = ".[" + str(modifications[site]) + "]" + modified_peptide
					else:
						raise click.ClickException("Error: Could not annotate N-terminus from peptide %s with delta mass %s." % (peptide['peptide_sequence'], nterm_modification))
				else:
					modified_peptide = ".(UniMod:" + str(record_id_nterm) + ")" + modified_peptide

			if cterm_modification is not "":
				record_id_cterm = um.get_id("C-term", 'Any C-term', cterm_modification)
				if record_id_cterm == -1:
					record_id_cterm = um.get_id("C-term", 'Protein C-term', cterm_modification)

				if record_id_cterm == -1:
					if self.enable_unannotated:
						modified_peptide = modified_peptide + ".[" + str(modifications[site]) + "]"
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
				base_name = posixpath.basename(ntpath.basename(elem.attrib['base_name']))

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

								if "pep" not in scores:
									# If 2 search hits have the same rank only the first one has the analysis_result explicitly written out.
									scores["pep"] = prev_pep
                                                                        
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

        for p in peptides:

            #percolator probability
            scores["q_value"] = float(p.getHits()[0].getMetaValue('MS:1001491'))
            scores["pep"] = float(p.getHits()[0].getMetaValue('MS:1001491'))

            try:
                parsed_peptides.append({**{'run_id': self.base_name,
                                           'scan_id': int(str(p.getMetaValue("spectrum_reference")).split('scan=')[-1].strip("'")),
                                           'hit_rank': int(p.getHits()[0].getRank()),
                                           'massdiff': float(0),
                                           'precursor_charge': int(p.getHits()[0].getCharge()),
                                           'retention_time': float(p.getRT()),
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
                                           'scan_id': int(str(p.getMetaValue("spectrum_reference")).split('scan=')[-1].strip("'")),
                                           'hit_rank': int(p.getHits()[0].getRank()),
                                           'massdiff': float(0),
                                           'precursor_charge': int(p.getHits()[0].getCharge()),
                                           'retention_time': float(p.getRT()),
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
		ptms_site = self.ptms[site]
		search_multiple_positions = isinstance(position, (list, tuple))

		kvs = itertools.chain.from_iterable((((k, p), v) for k,v in ptms_site[p].items()) for p in position) \
			if search_multiple_positions else \
			ptms_site[position].items()
		for key, value in kvs:
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

def read_mzml_or_mzxml_impl(path, psms, theoretical, max_delta_ppm, filetype):
	assert filetype in ('mzml', 'mzxml')
	fh = po.MzMLFile() if filetype=='mzml' else po.MzXMLFile()
	fh.setLogType(po.LogType.CMD)
	input_map = po.MSExperiment()
	fh.load(path, input_map)

	peaks_list = []
	for ix, psm in psms.iterrows():
		scan_id = psm['scan_id']
		ionseries = theoretical[psm['modified_peptide']][psm['precursor_charge']]

		spectrum = input_map.getSpectrum(scan_id - 1)

		fragments = []
		product_mzs = []
		intensities = []
		for peak in spectrum:
			fragment, product_mz = annotate_mass(peak.getMZ(), ionseries, max_delta_ppm)
			if fragment is not None:
				fragments.append(fragment)
				product_mzs.append(product_mz)
				intensities.append(peak.getIntensity())

		peaks = pd.DataFrame({'fragment': fragments, 'product_mz': product_mzs, 'intensity': intensities})
		peaks['scan_id'] = scan_id
		peaks['precursor_mz'] = po.AASequence.fromString(po.String(psm['modified_peptide'])).getMonoWeight(po.Residue.ResidueType.Full, psm['precursor_charge']) / psm['precursor_charge'];
		peaks['modified_peptide'] = psm['modified_peptide']
		peaks['precursor_charge'] = psm['precursor_charge']

		# Baseline normalization to highest annotated peak
		max_intensity = np.max(peaks['intensity'])
		if max_intensity > 0:
			peaks['intensity'] = peaks['intensity'] * (10000 / max_intensity)

		peaks_list.append(peaks)

	if len(peaks_list) > 0:
		transitions = pd.concat(peaks_list)
		# Multiple peaks might be identically annotated, only use most intense
		transitions = transitions.groupby(['scan_id','modified_peptide','precursor_charge','precursor_mz','fragment','product_mz'])['intensity'].max().reset_index()
	else:
		transitions = pd.DataFrame({'scan_id': [], 'modified_peptide': [], 'precursor_charge': [], 'precursor_mz': [], 'fragment': [], 'product_mz': [], 'intensity': []})
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
		ionseries = theoretical[modified_peptide][precursor_charge]

		mz_intensity_array = tims_data[scan_id]

		fragments = []
		product_mzs = []
		intensities = []
		for mz, intensity in mz_intensity_array:
			fragment, product_mz = annotate_mass(mz, ionseries, max_delta_ppm)
			if fragment is not None:
				fragments.append(fragment)
				product_mzs.append(product_mz)
				intensities.append(intensity)

		peaks = pd.DataFrame({'fragment': fragments, 'product_mz': product_mzs, 'intensity': intensities})
		peaks['scan_id'] = scan_id
		peaks['precursor_mz'] = po.AASequence.fromString(po.String(modified_peptide)).getMonoWeight(po.Residue.ResidueType.Full, precursor_charge) / precursor_charge;
		peaks['modified_peptide'] = modified_peptide
		peaks['precursor_charge'] = precursor_charge

		# Baseline normalization to highest annotated peak
		peaks['intensity'] = peaks['intensity'] * (10000 / np.max(peaks['intensity']))

		peaks_list.append(peaks)

	if len(peaks_list) > 0:
		transitions = pd.concat(peaks_list)
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


def generate_ionseries(peptide_sequence, precursor_charge, fragment_charges=[1,2,3,4], fragment_types=['b','y'], enable_specific_losses = False, enable_unspecific_losses = False):
	peptide = po.AASequence.fromString(po.String(peptide_sequence))
	sequence = peptide.toUnmodifiedString()

	unspecific_losses = ["H2O1","H3N1","C1H2N2","C1H2N1O1"]

	fragments = {}

	for fragment_type in fragment_types:
		for fragment_charge in fragment_charges:
			if fragment_charge <= precursor_charge:
				for fragment_ordinal in range(1,len(sequence)):
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

					# Standard fragment ions
					fragments[fragment_type + str(fragment_ordinal) + "^" + str(fragment_charge)] = mass

					# Losses
					if enable_specific_losses or enable_unspecific_losses:
						for lossfragment_ordinal in range(1,ion.size()):
							if (ion.getResidue(lossfragment_ordinal).hasNeutralLoss()):
								losses = ion.getResidue(lossfragment_ordinal).getLossFormulas()
								for loss in losses:
									loss_type = loss.toString()

									if (enable_specific_losses and loss_type not in unspecific_losses) or (enable_unspecific_losses and loss_type in unspecific_losses):
										fragments[fragment_type + str(fragment_ordinal) + "-" + loss_type + "^" + str(fragment_charge)] = mass - (loss.getMonoWeight() / fragment_charge)

	return list(fragments.keys()), np.fromiter(fragments.values(), np.float, len(fragments))

def conversion(pepxmlfile, spectralfile, unimodfile, exclude_range, max_delta_unimod, max_delta_ppm, enable_unannotated, enable_massdiff, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses):
	# Parse basename
	base_name = basename_spectralfile(spectralfile)
	click.echo("Info: Parsing run %s." % base_name)

	# Initialize UniMod
	um = unimod(unimodfile, max_delta_unimod)

	# Parse pepXML or idXML
	if pepxmlfile.casefold().endswith(('.pepxml', '.pep.xml')):
		click.echo("Info: Parsing pepXML.")
		px = pepxml(pepxmlfile, um, base_name, exclude_range, enable_unannotated, enable_massdiff)
	elif pepxmlfile.lower().endswith('idxml'):
		click.echo("Info: Parsing idXML.")
		px = idxml(pepxmlfile, base_name)
	else:
		click.echo('unknown format of pepxml identification file')

	# Continue if any PSMS are present
	psms = px.get()

	if psms.shape[0] > 0:
		run_id = basename_spectralfile(spectralfile)
		rank = re.compile(r'_rank([0-9]+)\.').search(pathlib.Path(pepxmlfile).name)
		rank_str = '' if rank is None else '_rank' + rank.group(1)
		psms['group_id'] = psms['run_id'] + "_" + psms['scan_id'].astype(str) + rank_str

		# Generate theoretical spectra
		click.echo("Info: Generate theoretical spectra.")
		theoretical = {}
		for ix, peptide in psms[['modified_peptide','precursor_charge']].drop_duplicates().iterrows():
			if peptide['modified_peptide'] not in theoretical.keys():
				theoretical[peptide['modified_peptide']] = {}

			theoretical[peptide['modified_peptide']][peptide['precursor_charge']] = generate_ionseries(peptide['modified_peptide'], peptide['precursor_charge'], fragment_charges, fragment_types, enable_specific_losses, enable_unspecific_losses)

		# Generate spectrum dataframe
		click.echo("Info: Processing spectra from file %s." % spectralfile)
		if spectralfile.lower().endswith(".mzxml"):
			peaks = read_mzml_or_mzxml_impl(spectralfile, psms[['scan_id','modified_peptide','precursor_charge']], theoretical, max_delta_ppm, 'mzxml')
		elif spectralfile.casefold().endswith(".mzml"):
			peaks = read_mzml_or_mzxml_impl(spectralfile, psms[['scan_id', 'modified_peptide', 'precursor_charge']], theoretical, max_delta_ppm, 'mzml')
		elif spectralfile.lower().endswith(".mgf"):
			peaks = read_mgf(spectralfile, psms[['scan_id', 'modified_peptide', 'precursor_charge']], theoretical, max_delta_ppm)

		# Round floating numbers
		peaks = peaks.round(6)

		return psms, peaks
	else:
		return pd.DataFrame({'run_id': [], 'scan_id': [], 'hit_rank': [], 'massdiff': [], 'precursor_charge': [], 'retention_time': [], 'ion_mobility': [], 'peptide_sequence': [], 'modifications': [], 'nterm_modification': [], 'cterm_modification': [], 'protein_id': [], 'gene_id': [], 'num_tot_proteins': [], 'decoy': []}), pd.DataFrame({'scan_id': [], 'modified_peptide': [], 'precursor_charge': [], 'precursor_mz': [], 'fragment': [], 'product_mz': [], 'intensity': []})

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
