import os
import subprocess
import shutil
import sys

import pandas as pd
import re

DATA_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _run_cmdline(cmdline):
    try:
        out = subprocess.check_output(cmdline, shell=True, stderr=subprocess.STDOUT)
        return out.decode(errors="replace")
    except subprocess.CalledProcessError as error:
        out = error.output.decode() if error.output else ""
        print(out, file=sys.stderr)
        raise


def _run_convertsage(temp_folder, regtest):
    os.chdir(temp_folder)

    # Copy test files to temp directory
    shutil.copy(os.path.join(DATA_FOLDER, "results.sage.tsv"), temp_folder)
    shutil.copy(os.path.join(DATA_FOLDER, "matched_fragments.sage.tsv"), temp_folder)

    cmdline = (
        "easypqp convertsage --sage_psm results.sage.tsv "
        "--sage_fragments matched_fragments.sage.tsv"
    )

    out = _run_cmdline(cmdline)
    # Strip leading timestamps of the form 'YYYY-MM-DD HH:MM:SS - ' and
    # filter out pyopenms environment warnings which are non-deterministic
    cleaned_lines = []
    for line in out.splitlines():
        line = re.sub(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} - ", "", line)
        # Remove pyopenms/openms data-path warnings which leak local site-packages paths
        if re.search(r"pyopenms", line, flags=re.IGNORECASE) or re.search(
            r"OPENMS_DATA_PATH", line
        ):
            continue
        # Remove pandas/bottleneck version warnings which are environment-specific
        if re.search(
            r"Pandas requires version", line, flags=re.IGNORECASE
        ) or re.search(r"bottleneck", line, flags=re.IGNORECASE):
            continue
        # Remove any raw warnings.warn(...) lines that may appear in some envs
        if "warnings.warn" in line:
            continue
        cleaned_lines.append(line)
    cleaned = "\n".join(cleaned_lines)
    print(cleaned, file=regtest)

    # Expect output files for run 'LQSRPAAPPAPGPGQLTLR'
    run_stem = "LQSRPAAPPAPGPGQLTLR"
    psmpkl = f"{run_stem}.psmpkl"
    peakpkl = f"{run_stem}.peakpkl"

    assert os.path.exists(psmpkl), f"Missing expected output {psmpkl}"
    assert os.path.exists(peakpkl), f"Missing expected output {peakpkl}"

    # Verify pickles load and contain expected columns
    psms = pd.read_pickle(psmpkl)
    peaks = pd.read_pickle(peakpkl)

    assert not psms.empty, "psmpkl is empty"
    assert not peaks.empty, "peakpkl is empty"

    # Basic schema checks
    assert "run_id" in psms.columns
    assert "scan_id" in psms.columns
    assert "run_id" in peaks.columns
    assert "product_mz" in peaks.columns or "fragment" in peaks.columns


def test_convertsage(tmpdir, regtest):
    _run_convertsage(tmpdir.strpath, regtest)
