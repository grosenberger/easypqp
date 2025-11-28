import os
import re
import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, List

from .util import timestamped_echo
from .convert import unimod as UniModHelper


def _basename_wo_ext(p: str) -> str:
    """
    Return the basename without extensions, handling common compression/archive
    extensions so that e.g. 'file.mzML.gz' -> 'file'.

    Strategy:
    - take the basename
    - if the final extension is a known compression/archive suffix ('.gz', '.bz2',
      '.zst', '.tgz', etc.) remove it
    - then remove one remaining extension (the common data file extension)
    - return the resulting stem
    """
    name = os.path.basename(p or "")
    comp_suffixes = {
        ".gz",
        ".zst",
        ".tar",
    }

    root, ext = os.path.splitext(name)
    if ext and ext.lower() in comp_suffixes:
        name = root

    stem, _ = os.path.splitext(name)
    return stem


def _get_first_existing(df: pd.DataFrame, cols: List[str], cast=None, default=None):
    for c in cols:
        if c in df.columns:
            return df[c] if cast is None else pd.to_numeric(df[c], errors="coerce")
    if default is None:
        return None
    return pd.Series([default] * len(df))


def _read_table(path: str) -> pd.DataFrame:
    """Read a TSV or Parquet file into a DataFrame."""
    p = (path or "").lower()
    if p.endswith((".parquet", ".pq")):
        try:
            return pd.read_parquet(path)
        except Exception as e:
            raise RuntimeError(f"Failed to read Parquet file: {path}\n{e}")
    if p.endswith((".tsv")):
        try:
            return pd.read_csv(path, sep="\t", dtype=str)
        except Exception as e:
            raise RuntimeError(f"Failed to read TSV file: {path}\n{e}")


class SagePSMParser:
    """
    Parse results.sage.tsv to EasyPQP PSM schema (subset used by library.generate)

    Output columns:
      run_id, scan_id, hit_rank, massdiff, precursor_charge, retention_time,
      ion_mobility, peptide_sequence, protein_id, gene_id, num_tot_proteins,
      decoy, pep, modified_peptide, group_id, precursor_mz
    """

    PROTON = 1.0072764
    NEUTRON = 1.00335
    # Sage bracket delta pattern: A[+15.9949], C[-0.9840], etc.
    BRACKET_RE = re.compile(r"([A-Z])\[(?P<delta>[+-]?\d+(?:\.\d+)?)\]")
    # Uniprot token pattern: db|ACCESSION|ENTRY_NAME  (e.g., sp|P01903|DRA_HUMAN)
    _ACC_ENTRY_RE = re.compile(r"^[A-Za-z]{2}\|(?P<acc>[^|]+)\|(?P<entry>[^|]+)$")
    # Common decoy prefixes occasionally carried into protein tokens (we still rely on label for decoy)
    _DECOY_PREFIX_RE = re.compile(r"^(?:decoy_|rev_)+", flags=re.IGNORECASE)

    def __init__(
        self,
        results_tsv: str,
        unimod_xml: Optional[str],
        max_delta_unimod: float = 0.02,
        mz_precision_digits: int = 6,
    ):
        self.results_tsv = results_tsv
        self.um = UniModHelper(unimod_xml, max_delta_unimod) if unimod_xml else None
        self.max_delta_unimod = max_delta_unimod
        self.mz_precision_digits = mz_precision_digits

    @staticmethod
    def _uniq_preserve(seq):
        """De-duplicate while preserving order."""
        seen = set()
        out = []
        for x in seq:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    def _clean_token(self, tok: str) -> str:
        """Strip decoy prefixes and whitespace from an individual protein token."""
        tok = (tok or "").strip()
        return self._DECOY_PREFIX_RE.sub("", tok)

    def _parse_protein_token(self, tok: str) -> Tuple[str, str]:
        """
        Extract (accession, entry_name) from a single token.
        Falls back gracefully if format isn't db|ACC|ENTRY.
        """
        t = self._clean_token(tok)
        m = self._ACC_ENTRY_RE.match(t)
        if m:
            return m.group("acc"), m.group("entry")
        # Fallbacks:
        if "|" in t:
            parts = t.split("|")
            if len(parts) >= 3:
                return parts[1] or "", parts[2] or ""
            # unknown pipe-y format: best-effort
            return parts[-1] or "", ""
        # No pipes at all: treat token as accession-only
        return t, ""

    def _split_accessions_and_entries(
        self, proteins: pd.Series
    ) -> Tuple[pd.Series, pd.Series, pd.Series]:
        """
        Vectorized split of Sage protein strings into:
          - accessions (semicolon-joined)
          - entry names (semicolon-joined)
          - count of unique accessions (for num_tot_proteins)
        """
        acc_list = []
        entry_list = []
        counts = []
        for s in proteins.astype(str):
            if not s or s == "nan":
                accs, entries = [], []
            else:
                toks = [t for t in s.split(";") if t.strip()]
                pairs = [self._parse_protein_token(t) for t in toks]
                accs = self._uniq_preserve([a for a, _ in pairs if a])
                entries = self._uniq_preserve([e for _, e in pairs if e])

            acc_list.append(";".join(accs))
            entry_list.append(";".join(entries))
            counts.append(len(accs))

        return pd.Series(acc_list), pd.Series(entry_list), pd.Series(counts)

    def _annotate_unimod(self, pep: str) -> str:
        """
        Convert Sage bracket deltas (e.g., M[+15.9949]) to (UniMod:<ID>).
        Tries position-specific contexts (N-term / C-term) before 'Anywhere'.
        Falls back to leaving the numeric delta if nothing matches.
        """
        if self.um is None or "[" not in pep:
            return pep

        # 1) get clean sequence and site->delta map from Sage string
        seq = re.sub(r"\[[-+0-9.]+\]", "", pep)
        site2delta: Dict[int, float] = {}
        site = 0
        i = 0
        while i < len(pep):
            ch = pep[i]
            if ch.isalpha():
                site += 1
                i += 1
                if i < len(pep) and pep[i] == "[":
                    j = pep.find("]", i + 1)
                    site2delta[site] = float(pep[i + 1 : j])
                    i = j + 1
            else:
                i += 1

        # 2) position preference helper
        def positions_for_site(idx: int, length: int):
            if idx == 1:
                # try N-terminus flavors first, then Anywhere
                return ["Any N-term", "Protein N-term", "Anywhere"]
            if idx == length:
                # try C-terminus flavors first, then Anywhere
                return ["Any C-term", "Protein C-term", "Anywhere"]
            return ["Anywhere"]

        # 3) very small fallback table for the most common N-term losses
        #    (used only if UniMod lookup fails)
        def fallback_unimod(aa: str, idx: int, delta: float, tol=0.02) -> int:
            if idx == 1 and aa == "Q" and abs(delta - (-17.026549)) <= tol:
                return 28  # Gln->pyro-Glu (N-term)
            if idx == 1 and aa == "E" and abs(delta - (-18.010565)) <= tol:
                return 27  # Glu->pyro-Glu (N-term)
            return -1

        # 4) build output by injecting (UniMod:<id>) after the modified residue
        out = list(seq)
        L = len(seq)
        for idx in sorted(site2delta.keys(), reverse=True):
            delta = site2delta[idx]
            aa = seq[idx - 1]
            rec_id = -1

            # Try position-specific contexts first
            for pos in positions_for_site(idx, L):
                rid = self.um.get_id(aa, pos, delta)
                if isinstance(rid, tuple):
                    rid = rid[0]
                if rid != -1:
                    rec_id = rid
                    break

            # Fallback: known N-term conversions (pyro-Glu/Q,E)
            if rec_id == -1:
                rec_id = fallback_unimod(aa, idx, delta, self.max_delta_unimod)

            insert = f"(UniMod:{rec_id})" if rec_id != -1 else f"[{delta:+.6f}]"
            out.insert(idx, insert)

        return "".join(out)

    def parse(self) -> pd.DataFrame:
        df = _read_table(self.results_tsv).fillna("")

        filename = _get_first_existing(
            df, ["filename", "file", "rawfile", "raw_file", "source_file"]
        )
        if filename is None:
            raise ValueError("results.sage.tsv is missing a filename/raw file column.")
        run_id = filename.astype(str).apply(_basename_wo_ext)

        scan_id = (
            _get_first_existing(
                df,
                ["scannr", "scan", "scan_id", "spectrum_index"],
                cast=float,
                default=np.nan,
            )
            .fillna(1)
            .astype(int)
        )
        hit_rank = (
            _get_first_existing(df, ["rank", "hit_rank"], cast=float, default=1)
            .fillna(1)
            .astype(int)
        )
        z = (
            _get_first_existing(
                df, ["precursor_charge", "charge", "z"], cast=float, default=2
            )
            .fillna(2)
            .astype(int)
        )

        rt = _get_first_existing(
            df,
            ["rt", "retention_time", "retention", "retention_time_sec"],
            cast=float,
            default=np.nan,
        )
        im = _get_first_existing(
            df, ["ion_mobility", "mobility", "ccs", "k0"], cast=float, default=np.nan
        )
        # If im is all 0s, set to NaN
        if im.eq(0).all():
            im = pd.Series([np.nan] * len(df))

        pep_seq = df["peptide"].astype(str)
        proteins_raw = _get_first_existing(df, ["proteins", "protein", "protein_id"])
        proteins_raw = (
            proteins_raw.astype(str)
            if proteins_raw is not None
            else pd.Series([""] * len(df))
        )
        protein_ids, gene_ids, num_prot = self._split_accessions_and_entries(
            proteins_raw
        )

        if "label" in df.columns:
            # decoy detection from label
            # Sage TSV: label == -1 (decoy), +1 (target)
            label_series = pd.to_numeric(df["label"], errors="coerce")
            decoy = label_series == -1
        elif "is_decoy" in df.columns:
            # The parquet format uses a boolean is_decoy column
            decoy = df["is_decoy"]

        # spectrum-level q-value, peptide-level q-value and protein-level q-value
        pep = (
            pd.to_numeric(df["posterior_error"], errors="coerce")
            if "posterior_error" in df.columns
            else pd.Series([np.nan] * len(df))
        )
        spectrum_q = (
            pd.to_numeric(df["spectrum_q"], errors="coerce")
            if "spectrum_q" in df.columns
            else pd.Series([np.nan] * len(df))
        )
        peptide_q = (
            pd.to_numeric(df["peptide_q"], errors="coerce")
            if "peptide_q" in df.columns
            else pd.Series([np.nan] * len(df))
        )
        protein_q = (
            pd.to_numeric(df["protein_q"], errors="coerce")
            if "protein_q" in df.columns
            else pd.Series([np.nan] * len(df))
        )

        # compute precursor m/z from neurtal mass using the theoretical calculated mass of the peptide.
        calcmass = _get_first_existing(df, ["calcmass"], cast=float, default=np.nan)
        prec_mz = pd.Series(np.nan, index=df.index, dtype=float)
        mask_calc = calcmass.notna() & (z > 0)
        prec_mz.loc[mask_calc] = (calcmass[mask_calc] + z[mask_calc] * self.PROTON) / z[
            mask_calc
        ]

        ## If we wanted to compute from experimental mass instead:
        # expmass  = _get_first_existing(df, ['expmass'],  cast=float, default=np.nan)
        # iso_err = _get_first_existing(
        #     df, ["isotope_error", "isotope"], cast=float, default=0.0
        # ).fillna(0.0)
        # mask_exp = prec_mz.isna() & expmass.notna() & (z > 0)
        # mz_exp = (expmass[mask_exp] + z[mask_exp] * PROTON) / z[mask_exp]
        # prec_mz.loc[mask_exp] = mz_exp - (iso_err[mask_exp] * NEUTRON) / z[mask_exp]

        ## set precision
        prec_mz = prec_mz.round(self.mz_precision_digits)

        # modified peptide
        modpep = pep_seq.apply(self._annotate_unimod)

        # group id (same style as convert paths)
        group_id = (
            run_id
            + "_"
            + scan_id.astype(str)
            + np.where(hit_rank > 1, "_rank" + hit_rank.astype(str), "")
        )

        out = pd.DataFrame(
            {
                "run_id": run_id,
                "scan_id": scan_id,
                "hit_rank": hit_rank,
                "massdiff": 0.0,
                "precursor_charge": z,
                "retention_time": rt,
                "ion_mobility": im,
                "peptide_sequence": pep_seq.str.replace(
                    r"\[[-+0-9.]+\]", "", regex=True
                ),
                "protein_id": protein_ids.fillna(""),
                "gene_id": gene_ids.fillna(""),
                "num_tot_proteins": num_prot.fillna(0).astype(int),
                "decoy": decoy.astype(bool),
                "modified_peptide": modpep,
                "group_id": group_id,
                "precursor_mz": prec_mz,
                "pep": pep,
                "q_value": spectrum_q,
                "peptide_q": peptide_q,
                "protein_q": protein_q,
            }
        )
        return out

    def parse_df(
        self, df: pd.DataFrame, psm_id_series: Optional[pd.Series] = None
    ) -> pd.DataFrame:
        """
        Parse a provided DataFrame slice (same logic as `parse` but works on an
        already-loaded DataFrame). This is useful for chunked/streaming flows.

        If `psm_id_series` is provided it will be attached to the returned
        DataFrame as a `psm_id` column (preserving positional alignment).
        """
        df = df.fillna("")

        filename = _get_first_existing(
            df, ["filename", "file", "rawfile", "raw_file", "source_file"]
        )
        if filename is None:
            raise ValueError("results.sage.tsv is missing a filename/raw file column.")
        run_id = filename.astype(str).apply(_basename_wo_ext)

        scan_id = (
            _get_first_existing(
                df,
                ["scannr", "scan", "scan_id", "spectrum_index"],
                cast=float,
                default=np.nan,
            )
            .fillna(1)
            .astype(int)
        )
        hit_rank = (
            _get_first_existing(df, ["rank", "hit_rank"], cast=float, default=1)
            .fillna(1)
            .astype(int)
        )
        z = (
            _get_first_existing(
                df, ["precursor_charge", "charge", "z"], cast=float, default=2
            )
            .fillna(2)
            .astype(int)
        )

        rt = _get_first_existing(
            df,
            ["rt", "retention_time", "retention", "retention_time_sec"],
            cast=float,
            default=np.nan,
        )
        im = _get_first_existing(
            df, ["ion_mobility", "mobility", "ccs", "k0"], cast=float, default=np.nan
        )
        # If im is all 0s, set to NaN
        if im.eq(0).all():
            im = pd.Series([np.nan] * len(df))

        pep_seq = df["peptide"].astype(str)
        proteins_raw = _get_first_existing(df, ["proteins", "protein", "protein_id"])
        proteins_raw = (
            proteins_raw.astype(str)
            if proteins_raw is not None
            else pd.Series([""] * len(df))
        )
        protein_ids, gene_ids, num_prot = self._split_accessions_and_entries(
            proteins_raw
        )

        if "label" in df.columns:
            label_series = pd.to_numeric(df["label"], errors="coerce")
            decoy = label_series == -1
        elif "is_decoy" in df.columns:
            decoy = df["is_decoy"]
        else:
            decoy = pd.Series([False] * len(df))

        pep = (
            pd.to_numeric(df["posterior_error"], errors="coerce")
            if "posterior_error" in df.columns
            else pd.Series([np.nan] * len(df))
        )
        spectrum_q = (
            pd.to_numeric(df["spectrum_q"], errors="coerce")
            if "spectrum_q" in df.columns
            else pd.Series([np.nan] * len(df))
        )
        peptide_q = (
            pd.to_numeric(df["peptide_q"], errors="coerce")
            if "peptide_q" in df.columns
            else pd.Series([np.nan] * len(df))
        )
        protein_q = (
            pd.to_numeric(df["protein_q"], errors="coerce")
            if "protein_q" in df.columns
            else pd.Series([np.nan] * len(df))
        )

        calcmass = _get_first_existing(df, ["calcmass"], cast=float, default=np.nan)
        prec_mz = pd.Series(np.nan, index=df.index, dtype=float)
        mask_calc = calcmass.notna() & (z > 0)
        prec_mz.loc[mask_calc] = (calcmass[mask_calc] + z[mask_calc] * self.PROTON) / z[
            mask_calc
        ]
        prec_mz = prec_mz.round(self.mz_precision_digits)

        modpep = pep_seq.apply(self._annotate_unimod)

        group_id = (
            run_id
            + "_"
            + scan_id.astype(str)
            + np.where(hit_rank > 1, "_rank" + hit_rank.astype(str), "")
        )

        out = pd.DataFrame(
            {
                "run_id": run_id,
                "scan_id": scan_id,
                "hit_rank": hit_rank,
                "massdiff": 0.0,
                "precursor_charge": z,
                "retention_time": rt,
                "ion_mobility": im,
                "peptide_sequence": pep_seq.str.replace(
                    r"\[[-+0-9.]+\]", "", regex=True
                ),
                "protein_id": protein_ids.fillna(""),
                "gene_id": gene_ids.fillna(""),
                "num_tot_proteins": num_prot.fillna(0).astype(int),
                "decoy": decoy.astype(bool),
                "modified_peptide": modpep,
                "group_id": group_id,
                "precursor_mz": prec_mz,
                "pep": pep,
                "q_value": spectrum_q,
                "peptide_q": peptide_q,
                "protein_q": protein_q,
            }
        )

        if psm_id_series is not None:
            out = out.reset_index(drop=True)
            out["psm_id"] = psm_id_series.astype(str).str.strip().reset_index(drop=True)

        return out


class SageFragmentParser:
    """
    Parse matched_fragments.sage.tsv to EasyPQP 'peaks' table:
      columns: scan_id, modified_peptide, precursor_charge, precursor_mz, fragment, product_mz, intensity
    """

    def __init__(self, frags_tsv: str, mz_precision_digits: int = 6):
        self.frags_tsv = frags_tsv
        self.mz_precision_digits = mz_precision_digits

    @staticmethod
    def _ann(ftype: str, ord_: int, z: int) -> str:
        return f"{ftype}{ord_}^{z}"

    def parse(self, psms_with_psmid: pd.DataFrame) -> pd.DataFrame:
        fr = _read_table(self.frags_tsv).fillna("")

        for c in [
            "psm_id",
            "fragment_ordinals",
            "fragment_charge",
            "fragment_mz_calculated",
            "fragment_mz_experimental",
            "fragment_intensity",
        ]:
            if c in fr.columns:
                fr[c] = pd.to_numeric(fr[c], errors="coerce")
        if "psm_id" not in fr.columns:
            raise ValueError(
                "matched_fragments.sage.tsv must contain a 'psm_id' column."
            )

        fr["psm_id"] = fr["psm_id"].astype(str).str.strip()

        fr["fragment"] = fr.apply(
            lambda r: self._ann(
                str(r["fragment_type"]),
                int(r["fragment_ordinals"]),
                int(r["fragment_charge"]),
            ),
            axis=1,
        )
        fr["product_mz"] = fr["fragment_mz_calculated"]

        # join to PSMs
        join_cols = [
            "psm_id",
            "scan_id",
            "modified_peptide",
            "precursor_mz",
            "precursor_charge",
            "run_id",
        ]
        j = fr.merge(psms_with_psmid[join_cols], on="psm_id", how="inner")

        peaks = j[
            [
                "run_id",
                "scan_id",
                "modified_peptide",
                "precursor_charge",
                "precursor_mz",
                "fragment",
                "product_mz",
                "fragment_intensity",
            ]
        ].copy()
        peaks.rename(columns={"fragment_intensity": "intensity"}, inplace=True)

        # per-PSM normalization to 10,000 (matches convert paths)
        peaks["intensity"] = peaks["intensity"].fillna(0.0)
        grp = peaks.groupby(
            ["run_id", "scan_id", "modified_peptide", "precursor_charge"], dropna=False
        )["intensity"]
        denom = grp.transform(lambda x: np.nanmax(x.values) if len(x) else np.nan)
        peaks["intensity"] = (peaks["intensity"] / denom) * 10000.0
        peaks["intensity"] = peaks["intensity"].fillna(0.0)

        # round and de-duplicate (keep most intense per exact fragment/product_mz)
        peaks["product_mz"] = peaks["product_mz"].round(self.mz_precision_digits)
        peaks["precursor_mz"] = peaks["precursor_mz"].round(self.mz_precision_digits)
        peaks["intensity"] = peaks["intensity"].round(self.mz_precision_digits)

        peaks = peaks.groupby(
            [
                "run_id",
                "scan_id",
                "modified_peptide",
                "precursor_charge",
                "precursor_mz",
                "fragment",
                "product_mz",
            ],
            as_index=False,
        )["intensity"].max()
        return peaks

    def parse_df(self, fr: pd.DataFrame, psms_with_psmid: pd.DataFrame) -> pd.DataFrame:
        """
        Parse a fragments DataFrame (filtered to relevant rows) and join to the
        provided PSM DataFrame. This mirrors `parse` but operates on in-memory
        DataFrames to support streaming.
        """
        fr = fr.fillna("")
        for c in [
            "psm_id",
            "fragment_ordinals",
            "fragment_charge",
            "fragment_mz_calculated",
            "fragment_mz_experimental",
            "fragment_intensity",
        ]:
            if c in fr.columns:
                fr[c] = pd.to_numeric(fr[c], errors="coerce")
        if "psm_id" not in fr.columns:
            raise ValueError(
                "matched_fragments.sage.tsv must contain a 'psm_id' column."
            )

        fr["psm_id"] = fr["psm_id"].astype(str).str.strip()

        fr["fragment"] = fr.apply(
            lambda r: self._ann(
                str(r["fragment_type"]),
                int(r["fragment_ordinals"]),
                int(r["fragment_charge"]),
            ),
            axis=1,
        )
        fr["product_mz"] = fr["fragment_mz_calculated"]

        join_cols = [
            "psm_id",
            "scan_id",
            "modified_peptide",
            "precursor_mz",
            "precursor_charge",
            "run_id",
        ]
        j = fr.merge(psms_with_psmid[join_cols], on="psm_id", how="inner")

        peaks = j[
            [
                "run_id",
                "scan_id",
                "modified_peptide",
                "precursor_charge",
                "precursor_mz",
                "fragment",
                "product_mz",
                "fragment_intensity",
            ]
        ].copy()
        peaks.rename(columns={"fragment_intensity": "intensity"}, inplace=True)

        peaks["intensity"] = peaks["intensity"].fillna(0.0)
        grp = peaks.groupby(
            ["run_id", "scan_id", "modified_peptide", "precursor_charge"], dropna=False
        )["intensity"]
        denom = grp.transform(lambda x: np.nanmax(x.values) if len(x) else np.nan)
        peaks["intensity"] = (peaks["intensity"] / denom) * 10000.0
        peaks["intensity"] = peaks["intensity"].fillna(0.0)

        peaks["product_mz"] = peaks["product_mz"].round(self.mz_precision_digits)
        peaks["precursor_mz"] = peaks["precursor_mz"].round(self.mz_precision_digits)
        peaks["intensity"] = peaks["intensity"].round(self.mz_precision_digits)

        peaks = peaks.groupby(
            [
                "run_id",
                "scan_id",
                "modified_peptide",
                "precursor_charge",
                "precursor_mz",
                "fragment",
                "product_mz",
            ],
            as_index=False,
        )["intensity"].max()
        return peaks


def convert_sage(
    results_tsv: str,
    fragments_tsv: str,
    unimod_xml: Optional[str],
    max_delta_unimod: float = 0.02,
    mz_precision_digits: int = 6,
    *,
    force_streaming: Optional[bool] = None,
    streaming_threshold_bytes: int = 1_000_000_000,
) -> Optional[List[str]]:
    """
    High-level conversion: Sage TSV/Parquet to EasyPQP PSM and peaks pickles written to disk.
    """
    # Auto-switch to streaming mode when inputs are very large, unless caller
    # explicitly requested non-streaming via force_streaming=False.
    try:
        if force_streaming is None:
            # determine combined size (fall back to streaming if either file is missing)
            try:
                rsize = os.path.getsize(results_tsv)
            except Exception:
                rsize = 0
            try:
                fsize = os.path.getsize(fragments_tsv)
            except Exception:
                fsize = 0
            use_stream = (rsize + fsize) >= streaming_threshold_bytes
        else:
            use_stream = bool(force_streaming)
    except Exception:
        use_stream = False

    if use_stream:
        timestamped_echo("Info: Using streaming Sage conversion for low memory usage")
        return convert_sage_streaming(
            results_tsv,
            fragments_tsv,
            unimod_xml,
            max_delta_unimod=max_delta_unimod,
            mz_precision_digits=mz_precision_digits,
        )

    # Read raw to extract psm_id for joining
    timestamped_echo("Info: Reading Sage PSMs")
    raw_res = _read_table(results_tsv)
    if "psm_id" not in raw_res.columns:
        raise ValueError(
            "results.sage.tsv must contain a 'psm_id' for joining with matched fragments."
        )

    raw_res["psm_id"] = raw_res["psm_id"].astype(str).str.strip()

    psms = SagePSMParser(
        results_tsv, unimod_xml, max_delta_unimod, mz_precision_digits
    ).parse()
    psms = raw_res[["psm_id"]].join(psms)

    if psms.empty:
        raise ValueError("No PSMs were parsed from the provided results.sage.tsv file.")

    timestamped_echo("Info: Reading Sage matched fragment peaks")
    peaks = SageFragmentParser(fragments_tsv, mz_precision_digits).parse(psms)

    if peaks.empty:
        raise ValueError(
            "No fragment peaks were parsed from the provided matched_fragments.sage.tsv file."
        )

    # Trim to minimal schema expected by library.generate
    keep = [
        "run_id",
        "scan_id",
        "hit_rank",
        "massdiff",
        "precursor_charge",
        "retention_time",
        "ion_mobility",
        "peptide_sequence",
        "protein_id",
        "gene_id",
        "num_tot_proteins",
        "decoy",
        "modified_peptide",
        "group_id",
        "pep",
        "q_value",
        "peptide_q",
        "protein_q",
    ]
    psms_export = psms[keep].copy()

    runs = sorted(psms_export["run_id"].dropna().unique().tolist())
    new_infiles = []
    for run in runs:
        psms_r = psms_export.loc[psms_export["run_id"] == run]
        peaks_r = (
            peaks.loc[peaks["run_id"] == run] if "run_id" in peaks.columns else peaks
        )

        if psms_r.empty or peaks_r.empty:
            timestamped_echo(
                f"Info: Skipping run {run}: psms={len(psms_r)}, peaks={len(peaks_r)}"
            )
            continue

        psmpkl = f"{run}.psmpkl"
        peakpkl = f"{run}.peakpkl"
        psms_r.to_pickle(psmpkl)
        peaks_r.to_pickle(peakpkl)
        timestamped_echo(f"Info: Wrote {psmpkl} and {peakpkl}")
        new_infiles.extend([psmpkl, peakpkl])

    if len(new_infiles) == 0:
        # click may not be available in all contexts; raise a generic error here
        raise RuntimeError("No non-empty runs detected after Sage conversion.")


def convert_sage_streaming(
    results_tsv: str,
    fragments_tsv: str,
    unimod_xml: Optional[str],
    max_delta_unimod: float = 0.02,
    mz_precision_digits: int = 6,
    chunksize: int = 100_000,
    tmpdir: Optional[str] = None,
) -> List[str]:
    """
    Memory-efficient streaming conversion that processes one run at a time.

    Returns a list of produced filenames (psmpkl and peakpkl pairs).
    """
    import tempfile

    # determine filename column quickly from header
    header = pd.read_csv(results_tsv, sep="\t", nrows=0)
    cols = header.columns.tolist()
    filename_col = None
    for c in ["filename", "file", "source_file"]:
        if c in cols:
            filename_col = c
            break
    if filename_col is None:
        raise ValueError("results.sage.tsv is missing a filename column.")
    if "psm_id" not in cols:
        raise ValueError(
            "results.sage.tsv must contain a 'psm_id' for joining with matched fragments."
        )

    runs = set()
    for chunk in pd.read_csv(
        results_tsv, sep="\t", usecols=[filename_col], dtype=str, chunksize=chunksize
    ):
        runs.update(
            chunk[filename_col]
            .fillna("")
            .astype(str)
            .apply(_basename_wo_ext)
            .unique()
            .tolist()
        )

    runs = sorted([r for r in runs if r])
    if not runs:
        raise RuntimeError("No runs discovered in results file")

    timestamped_echo(f"Info: Discovered {len(runs)} runs")

    outfiles = []
    tmpdir = tmpdir or tempfile.mkdtemp(prefix="easypqp_sage_")

    # Process one run at a time: stream results, parse chunks for this run only,
    # then stream fragments and filter psm_ids for this run only.
    for run in runs:
        timestamped_echo(f"Info: Processing run {run}")

        # gather parsed PSMs for this run in small chunks
        psm_chunks = []
        for chunk in pd.read_csv(results_tsv, sep="\t", dtype=str, chunksize=chunksize):
            # compute run ids for this chunk
            filename = _get_first_existing(
                chunk, ["filename", "file", "rawfile", "raw_file", "source_file"]
            )
            if filename is None:
                filename = pd.Series([""] * len(chunk))
            chunk_run = filename.astype(str).apply(_basename_wo_ext)
            mask = chunk_run == run
            if not mask.any():
                continue
            filtered = chunk.loc[mask].copy()
            psm_id_series = (
                filtered["psm_id"].astype(str).str.strip()
                if "psm_id" in filtered.columns
                else None
            )
            parsed = SagePSMParser(
                results_tsv, unimod_xml, max_delta_unimod, mz_precision_digits
            ).parse_df(filtered, psm_id_series=psm_id_series)
            psm_chunks.append(parsed)

        if not psm_chunks:
            timestamped_echo(f"Info: Skipping run {run}: no PSMs")
            continue
        psms_run = pd.concat(psm_chunks, ignore_index=True)

        # Now stream fragments and collect those matching psms_run["psm_id"]
        psm_id_set = set(psms_run["psm_id"].astype(str).tolist())
        frag_chunks = []
        for fr_chunk in pd.read_csv(
            fragments_tsv, sep="\t", dtype=str, chunksize=chunksize
        ):
            if "psm_id" not in fr_chunk.columns:
                continue
            maskf = fr_chunk["psm_id"].astype(str).str.strip().isin(psm_id_set)
            if not maskf.any():
                continue
            frag_chunks.append(fr_chunk.loc[maskf].copy())

        if not frag_chunks:
            timestamped_echo(f"Info: Skipping run {run}: no fragment peaks")
            continue
        fr_run = pd.concat(frag_chunks, ignore_index=True)

        # parse fragment chunks joined to psms_run
        peaks_run = SageFragmentParser(fragments_tsv, mz_precision_digits).parse_df(
            fr_run, psms_run
        )

        # Trim to expected schema and write pickles
        keep = [
            "run_id",
            "scan_id",
            "hit_rank",
            "massdiff",
            "precursor_charge",
            "retention_time",
            "ion_mobility",
            "peptide_sequence",
            "protein_id",
            "gene_id",
            "num_tot_proteins",
            "decoy",
            "modified_peptide",
            "group_id",
            "pep",
            "q_value",
            "peptide_q",
            "protein_q",
        ]
        psms_export = psms_run[keep].copy()

        if psms_export.empty or peaks_run.empty:
            timestamped_echo(
                f"Info: Skipping run {run}: psms={len(psms_export)}, peaks={len(peaks_run)}"
            )
            continue

        psmpkl = f"{run}.psmpkl"
        peakpkl = f"{run}.peakpkl"
        psms_export.to_pickle(psmpkl)
        peaks_run.to_pickle(peakpkl)
        timestamped_echo(f"Info: Wrote {psmpkl} and {peakpkl}")
        outfiles.extend([psmpkl, peakpkl])

    if not outfiles:
        raise RuntimeError(
            "No non-empty runs detected after Sage streaming conversion."
        )

    return outfiles
