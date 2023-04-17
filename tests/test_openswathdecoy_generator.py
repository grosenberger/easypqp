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

def _run_openswath_decoy_generator(regtest, temp_folder):
    os.chdir(temp_folder)
    data_path = os.path.join(DATA_FOLDER, "library_targets.pqp")
    shutil.copy(data_path, temp_folder)
    cmdline = "easypqp openswath-decoy-generator  --in library_targets.pqp --out library.pqp --method pseudo-reverse"

    stdout = _run_cmdline(cmdline)

    conn = sqlite3.connect("library.pqp")
    protein_table = pd.read_sql_query("SELECT * FROM PROTEIN", conn)
    peptide_table = pd.read_sql_query("SELECT * FROM PEPTIDE", conn)
    precursor_table = pd.read_sql_query("SELECT * FROM PRECURSOR", conn)
    transition_table = pd.read_sql_query("SELECT * FROM TRANSITION", conn)
    conn.close()

    print(protein_table.sort_values("ID"),file=regtest)
    print(peptide_table.sort_values("ID"),file=regtest)
    print(precursor_table.sort_values("ID"),file=regtest)
    print(transition_table.sort_values("ID"),file=regtest)

def test_openswath_decoy_generator(tmpdir, regtest):
    _run_openswath_decoy_generator(regtest, tmpdir.strpath)