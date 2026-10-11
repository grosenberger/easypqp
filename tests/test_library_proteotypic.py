import os
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from easypqp.library import HAS_PYPROPHET, process_psms

DATA_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def run_library(tmpdir, outfile, *options):
    os.chdir(tmpdir.strpath)
    psmpkl = "HepG2_rep1_small.psmpkl"
    peakpkl = "HepG2_rep1_small.peakpkl"
    shutil.copy(os.path.join(DATA_FOLDER, psmpkl), tmpdir.strpath)
    shutil.copy(os.path.join(DATA_FOLDER, peakpkl), tmpdir.strpath)

    subprocess.check_call([
        "easypqp", "library",
        "--out", outfile,
        "--nofdr",
        "--perform_rt_calibration", "False",
        "--diannpqp",
        *options,
        psmpkl, peakpkl,
    ])

    return pd.read_csv(outfile, sep="\t")


def test_library_no_proteotypic(tmpdir, regtest):
    proteotypic = run_library(tmpdir, "proteotypic.tsv")
    shared = run_library(tmpdir, "no_proteotypic.tsv", "--no-proteotypic")

    # The default keeps only single-protein peptides
    assert not proteotypic["ProteinId"].str.contains("[;,]").any()
    assert (proteotypic["Proteotypic"] == 1).all()

    # --no-proteotypic adds the shared peptides, with all proteins joined by ";"
    assert set(proteotypic["ModifiedPeptideSequence"]) < set(shared["ModifiedPeptideSequence"])
    assert not shared["ProteinId"].str.contains(",").any()
    assert (shared["Proteotypic"] == (~shared["ProteinId"].str.contains(";")).astype(int)).all()
    assert (shared["Proteotypic"] == 0).any()

    precursors = shared[
        ["ModifiedPeptideSequence", "PrecursorCharge", "ProteinId", "Proteotypic"]
    ].drop_duplicates()
    precursors.sort_values(["ModifiedPeptideSequence", "PrecursorCharge"]).to_csv(
        regtest, sep="\t", index=False
    )


@pytest.mark.skipif(not HAS_PYPROPHET, reason="pyprophet is not installed")
@pytest.mark.parametrize("proteotypic", [True, False])
def test_process_psms_shared_peptide_fdr(tmpdir, proteotypic):
    def peptides(prefix, n, pp, protein_prefix, decoy):
        return pd.DataFrame({
            "peptide_sequence": [f"{prefix}{i}" for i in range(n)],
            "protein_id": [f"{protein_prefix}{i}" for i in range(n)],
            "num_tot_proteins": 1,
            "decoy": decoy,
            "pp": pp,
        })

    shared = pd.DataFrame({
        "peptide_sequence": ["SHAREDPEPTIDE", "SHAREDPEPTIDE"],
        "protein_id": ["PB;PA", "PA,PB"],
        "num_tot_proteins": 2,
        "decoy": False,
        "pp": 0.999,
    })
    psms = pd.concat([
        peptides("UNIQUE", 50, 0.999, "P", False),
        peptides("NOISE", 100, np.linspace(0.0, 0.5, 100), "N", False),
        peptides("DECOY", 1000, np.linspace(0.0, 0.5, 1000), "DECOY_N", True),
        shared,
    ], ignore_index=True)
    psms["modified_peptide"] = psms["peptide_sequence"]
    psms["pep"] = 1 - psms["pp"]
    psms["run_id"] = "run1.mzML"

    result = process_psms(
        psms, None, None, 0.01, 0.01, 0.01, np.arange(0.1, 0.5, 0.05),
        os.path.join(tmpdir.strpath, "peptide.pdf"),
        os.path.join(tmpdir.strpath, "protein.pdf"),
        proteotypic, False,
    )

    assert set(result["peptide_sequence"]) >= {f"UNIQUE{i}" for i in range(50)}
    assert not result["decoy"].any()
    if proteotypic:
        assert "SHAREDPEPTIDE" not in set(result["peptide_sequence"])
    else:
        shared_result = result[result["peptide_sequence"] == "SHAREDPEPTIDE"]
        assert shared_result.shape[0] == 2
        assert set(shared_result["protein_id"]) == {"PA;PB"}
