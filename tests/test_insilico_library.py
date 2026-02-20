from __future__ import print_function

import os
import subprocess
import shutil
import sys

import pandas as pd
import pytest

pd.options.display.expand_frame_repr = False
pd.options.display.precision = 4
pd.options.display.max_columns = None

DATA_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# Check if insilico feature is available
try:
    from easypqp_rs import generate_insilico_library  # noqa: F401

    HAS_RUST_BACKEND = True
except ImportError:
    HAS_RUST_BACKEND = False


def _run_cmdline(cmdline):
    stdout = cmdline + "\n"
    try:
        stdout += str(
            subprocess.check_output(cmdline, shell=True, stderr=subprocess.STDOUT)
        )
    except subprocess.CalledProcessError as error:
        print(error, end="", file=sys.stderr)
        print(
            "Command output:",
            error.output.decode() if error.output else "No output",
            file=sys.stderr,
        )
        raise
    return stdout


def _run_insilico_library(regtest, temp_folder):
    os.chdir(temp_folder)
    config_path = os.path.join(DATA_FOLDER, "config.json")
    fasta_path = os.path.join(DATA_FOLDER, "Q99536.fasta")

    # Copy test files to temp directory
    shutil.copy(config_path, temp_folder)
    shutil.copy(fasta_path, temp_folder)

    # Update config to use local paths in temp folder
    import json

    with open("config.json", "r") as f:
        config = json.load(f)

    # Update paths to be relative to temp folder
    config["database"]["fasta"] = "Q99536.fasta"
    config["output_file"] = "easypqp_insilico_library.tsv"

    with open("config.json", "w") as f:
        json.dump(config, f, indent=2)

    cmdline = "easypqp insilico-library --config config.json"

    _run_cmdline(cmdline)

    # Read and verify the output TSV file
    output_file = "easypqp_insilico_library.tsv"
    assert os.path.exists(output_file), f"Output file {output_file} was not created"

    library_df = pd.read_csv(output_file, sep="\t")

    # Print basic statistics about the generated library
    print(f"Generated library contains {len(library_df)} transitions", file=regtest)

    # Use appropriate column for peptide count (prefer ModifiedPeptideSequence)
    peptide_col = (
        "ModifiedPeptideSequence"
        if "ModifiedPeptideSequence" in library_df.columns
        else "PeptideSequence"
    )

    # Compute unique precursors. Prefer TransitionGroupId if present; otherwise
    # group by PrecursorMz, PrecursorCharge and the peptide column (which may
    # include modification information). This avoids conflating different
    # peptidoforms that share the same base sequence.
    if "TransitionGroupId" in library_df.columns:
        num_precursors = library_df["TransitionGroupId"].nunique()
    else:
        num_precursors = library_df.groupby(
            ["PrecursorMz", "PrecursorCharge", peptide_col]
        ).ngroups

    print(f"Number of unique precursors: {num_precursors}", file=regtest)
    print(
        f"Number of unique peptides: {library_df[peptide_col].nunique()}", file=regtest
    )

    # Print column names
    print(f"\nColumns: {list(library_df.columns)}", file=regtest)

    # Round LibraryIntensity to make test more stable (DL predictions can vary slightly)
    # Keep only deterministic columns for display. Sort the entire library first
    # to ensure a stable, deterministic sample (taking head() after sorting can
    # otherwise pick different rows across runs). Use a stable sort (mergesort)
    # and include several tie-breaker columns to avoid non-determinism when many
    # fragments share the same precursor/product m/z values.
    preferred_sort_cols = [
        "PrecursorMz",
        "PrecursorCharge",
        "ProductMz",
        "ProductCharge",
        peptide_col,
        "FragmentType",
        "FragmentSeriesNumber",
        "Annotation",
        "TransitionId",
    ]
    sort_cols = [c for c in preferred_sort_cols if c in library_df.columns]
    if not sort_cols:
        # Fallback to a minimal deterministic sort
        sort_cols = ["PrecursorMz", "ProductMz"]
    display_df = (
        library_df.sort_values(sort_cols, kind="mergesort")
        .reset_index(drop=True)
        .head()
        .copy()
    )
    if "LibraryIntensity" in display_df.columns:
        display_df["LibraryIntensity"] = display_df["LibraryIntensity"].round(0)

    # Sort by ProductMz to ensure consistent ordering across runs
    display_df = display_df.sort_values("ProductMz").reset_index(drop=True)

    # Print a sample of the data (first few rows) - excluding non-deterministic columns
    print("\nFirst 5 transitions (deterministic columns only):", file=regtest)
    deterministic_cols = [
        "PrecursorMz",
        "ProductMz",
        "PrecursorCharge",
        "ProductCharge",
        "PeptideSequence",
        "ModifiedPeptideSequence",
        "ProteinId",
        "UniprotId",
        "GeneName",
        "FragmentType",
        "FragmentSeriesNumber",
        "Annotation",
        "TransitionGroupId",
        "TransitionId",
        "Decoy",
    ]
    available_cols = [col for col in deterministic_cols if col in display_df.columns]
    print(display_df[available_cols].to_string(), file=regtest)

    # Verify core columns exist (using actual column names from the output)
    core_columns = [
        "PrecursorMz",
        "ProductMz",
        "PrecursorCharge",
        "ProductCharge",
        "LibraryIntensity",
        "PeptideSequence",
        "ProteinId",
        "FragmentType",
        "FragmentSeriesNumber",
        "Annotation",
    ]

    missing_columns = [col for col in core_columns if col not in library_df.columns]
    if missing_columns:
        print(f"\nWarning: Missing core columns: {missing_columns}", file=regtest)

    # Print some statistics
    print("\nStatistics:", file=regtest)
    print(
        f"  Precursor charge range: {library_df['PrecursorCharge'].min()}-{library_df['PrecursorCharge'].max()}",
        file=regtest,
    )
    print(
        f"  Fragment types: {sorted(library_df['FragmentType'].unique())}", file=regtest
    )

    # Check for decoys using the Decoy column
    if "Decoy" in library_df.columns:
        print(f"  Contains decoys: {library_df['Decoy'].sum() > 0}", file=regtest)
        print(f"  Number of targets: {(library_df['Decoy'] == 0).sum()}", file=regtest)
        print(f"  Number of decoys: {(library_df['Decoy'] == 1).sum()}", file=regtest)
    else:
        print("  Decoy column not found in output", file=regtest)

    # Verify LibraryIntensity values are reasonable (not in regtest output due to variance)
    if "LibraryIntensity" in library_df.columns:
        intensity_stats = library_df["LibraryIntensity"].describe()
        # Only assert, don't print to regtest to avoid flakiness
        assert intensity_stats["min"] >= 0, "LibraryIntensity should be non-negative"
        assert intensity_stats["max"] <= 10001, (
            "LibraryIntensity should be normalized to ~10000"
        )


@pytest.mark.skipif(
    not HAS_RUST_BACKEND,
    reason="In-silico feature not installed (easypqp_rs package missing - reinstall easypqp)",
)
def test_insilico_library(tmpdir, regtest):
    _run_insilico_library(regtest, tmpdir.strpath)
