import numpy as np
import pandas as pd
import sys

from Bio import SeqIO

def is_proteotypic(protein_id):
	if ";" in protein_id:
		return False
	else:
		return True

def is_decoy(protein_id):
	decoy = True
	for protein in protein_id.split(";"):
		if "DECOY" not in protein:
			decoy = False
	return decoy

def parse_fasta(file):
	peptides = {}

	for record in SeqIO.parse(file, "fasta"):
		if str(record.seq) not in peptides.keys():
			peptides[str(record.seq)] = []
		peptides[str(record.seq)].append(record.id)

	peptide_sequences = []
	protein_ids = []
	for key, value in peptides.items():
		peptide_sequences.append(key)
		protein_ids.append(";".join(set(value)))

	peptide_index = pd.DataFrame({'protein_id': protein_ids, 'peptide_sequence': peptide_sequences})

	peptide_index['proteotypic'] = peptide_index['protein_id'].apply(is_proteotypic)
	peptide_index['decoy'] = peptide_index['protein_id'].apply(is_decoy)

	return peptide_index
