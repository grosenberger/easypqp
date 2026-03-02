import os
import subprocess
import shutil

import pandas as pd

DATA_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def test_library_diannpqp(tmpdir, regtest):
    os.chdir(tmpdir.strpath)
    psmpkl = "HepG2_rep1_small.psmpkl"
    peakpkl = "HepG2_rep1_small.peakpkl"
    shutil.copy(os.path.join(DATA_FOLDER, psmpkl), tmpdir.strpath)
    shutil.copy(os.path.join(DATA_FOLDER, peakpkl), tmpdir.strpath)

    outfile = "diannpqp_test_output.tsv"
    subprocess.check_call([
        "easypqp", "library",
        "--out", outfile,
        "--nofdr",
        "--perform_rt_calibration", "False",
        "--diannpqp",
        psmpkl, peakpkl,
    ])

    df = pd.read_csv(outfile, sep="\t")

    # Verify DIA-NN2 columns exist
    for col in ['FragmentType', 'FragmentSeriesNumber', 'FragmentLossType', 'FragmentCharge', 'Proteotypic']:
        assert col in df.columns, f"Missing column: {col}"

    # Verify Annotation column was removed
    assert 'Annotation' not in df.columns

    # Verify FragmentSeriesNumber is correctly parsed for multi-digit values
    assert (df['FragmentSeriesNumber'] >= 10).any(), \
        "Expected multi-digit FragmentSeriesNumber values (>=10) but found none"

    # Verify FragmentType values are valid
    assert df['FragmentType'].isin(['b', 'y']).all()

    # Verify FragmentCharge values are positive integers
    assert (df['FragmentCharge'] > 0).all()

    # Write sorted output for regression testing (TSV for deterministic format)
    df.sort_values(['PrecursorMz', 'ProductMz']).reset_index(drop=True).to_csv(regtest, sep="\t", index=False)
