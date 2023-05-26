from __future__ import print_function

import os
import subprocess
import shutil
import sys

import pandas as pd
import sqlite3

import pytest

pd.options.display.expand_frame_repr = False
pd.options.display.precision = 4
pd.options.display.max_columns = None

DATA_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

def _run_cmdline(cmdline):
    stdout = cmdline + "\n"
    try:
        stdout += str(subprocess.check_output(cmdline, shell=True,
                                          stderr=subprocess.STDOUT))
    except subprocess.CalledProcessError as error:
        print(error, end="", file=sys.stderr)
        raise
    return stdout

def _run_targetedfileconverter(regtest, temp_folder, infile, outfile):
    os.chdir(temp_folder)
    data_path = os.path.join(DATA_FOLDER, infile)
    shutil.copy(data_path, temp_folder)
    cmdline = f"easypqp targeted-file-converter  --in {infile} --out {outfile}"

    stdout = _run_cmdline(cmdline)

    if outfile.split(".")[1] == "pqp":
        conn = sqlite3.connect(outfile)
        protein_table = pd.read_sql_query("SELECT * FROM PROTEIN", conn)
        peptide_table = pd.read_sql_query("SELECT * FROM PEPTIDE", conn)
        precursor_table = pd.read_sql_query("SELECT * FROM PRECURSOR", conn)
        transition_table = pd.read_sql_query("SELECT * FROM TRANSITION", conn)
        conn.close()

        print(protein_table.sort_values("ID"),file=regtest)
        print(peptide_table.sort_values("ID"),file=regtest)
        print(precursor_table.sort_values("ID"),file=regtest)
        print(transition_table.sort_values("ID"),file=regtest)
    elif outfile.split(".")[1] == "tsv":
        print(pd.read_csv(outfile, sep="\t", nrows=100).sort_index(axis=1),file=regtest)


def test_targeted_file_converter_tsvtopqp(tmpdir, regtest):
    _run_targetedfileconverter(regtest, tmpdir.strpath, "test_transition_list.tsv", "test_transition_list.pqp")

def test_targeted_file_converter_pqptotsv(tmpdir, regtest):
    _run_targetedfileconverter(regtest, tmpdir.strpath, "test_transition_list.pqp", "test_transition_list.tsv")