import os
import re
import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, List

from .util import timestamped_echo
from .convert import unimod as UniModHelper


def _basename_wo_ext(p: str) -> str:
    return os.path.splitext(os.path.basename(p))[0]


def _get_first_existing(df: pd.DataFrame, cols: List[str], cast=None, default=None):
    for c in cols:
        if c in df.columns:
            return df[c] if cast is None else pd.to_numeric(df[c], errors='coerce')
    if default is None:
        return None
    return pd.Series([default] * len(df))


class SagePSMParser:
    """
    Parse results.sage.tsv to EasyPQP PSM schema (subset used by library.generate)

    Output columns:
      run_id, scan_id, hit_rank, massdiff, precursor_charge, retention_time,
      ion_mobility, peptide_sequence, protein_id, gene_id, num_tot_proteins,
      decoy, pep, modified_peptide, group_id, precursor_mz (helper for join)
    """
    # Sage bracket delta pattern: A[+15.9949], C[-0.9840], etc.
    BRACKET_RE = re.compile(r'([A-Z])\[(?P<delta>[+-]?\d+(?:\.\d+)?)\]')
    # Uniprot token pattern: db|ACCESSION|ENTRY_NAME  (e.g., sp|P01903|DRA_HUMAN)
    _ACC_ENTRY_RE = re.compile(r'^[A-Za-z]{2}\|(?P<acc>[^|]+)\|(?P<entry>[^|]+)$')
    # Common decoy prefixes occasionally carried into protein tokens (we still rely on label for decoy)
    _DECOY_PREFIX_RE = re.compile(r'^(?:decoy_|rev_)+', flags=re.IGNORECASE)


    def __init__(self, results_tsv: str, unimod_xml: Optional[str], max_delta_unimod: float = 0.02):
        self.results_tsv = results_tsv
        self.um = UniModHelper(unimod_xml, max_delta_unimod) if unimod_xml else None
        
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
            return m.group('acc'), m.group('entry')
        # Fallbacks:
        if '|' in t:
            parts = t.split('|')
            if len(parts) >= 3:
                return parts[1] or '', parts[2] or ''
            # unknown pipe-y format: best-effort
            return parts[-1] or '', ''
        # No pipes at all: treat token as accession-only
        return t, ''

    def _split_accessions_and_entries(self, proteins: pd.Series) -> Tuple[pd.Series, pd.Series, pd.Series]:
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
            if not s or s == 'nan':
                accs, entries = [], []
            else:
                toks = [t for t in s.split(';') if t.strip()]
                pairs = [self._parse_protein_token(t) for t in toks]
                accs   = self._uniq_preserve([a for a, _ in pairs if a])
                entries= self._uniq_preserve([e for _, e in pairs if e])

            acc_list.append(';'.join(accs))
            entry_list.append(';'.join(entries))
            counts.append(len(accs))

        return pd.Series(acc_list), pd.Series(entry_list), pd.Series(counts)

    def _annotate_unimod(self, pep: str) -> str:
        """
        Convert Sage bracket deltas (e.g., M[+15.9949]) to (UniMod:<ID>).
        Tries position-specific contexts (N-term / C-term) before 'Anywhere'.
        Falls back to leaving the numeric delta if nothing matches.
        """
        if self.um is None or '[' not in pep:
            return pep

        # 1) get clean sequence and site->delta map from Sage string
        seq = re.sub(r'\[[-+0-9.]+\]', '', pep)
        site2delta: Dict[int, float] = {}
        site = 0
        i = 0
        while i < len(pep):
            ch = pep[i]
            if ch.isalpha():
                site += 1
                i += 1
                if i < len(pep) and pep[i] == '[':
                    j = pep.find(']', i + 1)
                    site2delta[site] = float(pep[i + 1:j])
                    i = j + 1
            else:
                i += 1

        # 2) position preference helper
        def positions_for_site(idx: int, length: int):
            if idx == 1:
                # try N-terminus flavors first, then Anywhere
                return ['Any N-term', 'Protein N-term', 'Anywhere']
            if idx == length:
                # try C-terminus flavors first, then Anywhere
                return ['Any C-term', 'Protein C-term', 'Anywhere']
            return ['Anywhere']

        # 3) very small fallback table for the most common N-term losses
        #    (used only if UniMod lookup fails)
        def fallback_unimod(aa: str, idx: int, delta: float) -> int:
            tol = 0.02
            if idx == 1 and aa == 'Q' and abs(delta - (-17.026549)) <= tol:
                return 28  # Gln->pyro-Glu (N-term)
            if idx == 1 and aa == 'E' and abs(delta - (-18.010565)) <= tol:
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
                rec_id = fallback_unimod(aa, idx, delta)

            insert = f"(UniMod:{rec_id})" if rec_id != -1 else f"[{delta:+.6f}]"
            out.insert(idx, insert)

        return ''.join(out)


    def parse(self) -> pd.DataFrame:
        df = pd.read_csv(self.results_tsv, sep='\t', dtype=str).fillna('')

        # basics
        filename = _get_first_existing(df, ['filename', 'file', 'rawfile', 'raw_file', 'source_file'])
        if filename is None:
            raise ValueError("results.sage.tsv is missing a filename/raw file column.")
        run_id = filename.astype(str).apply(_basename_wo_ext)

        scan_id = _get_first_existing(df, ['scannr', 'scan', 'scan_id', 'spectrum_index'], cast=float, default=np.nan).fillna(1).astype(int)
        hit_rank = _get_first_existing(df, ['rank', 'hit_rank'], cast=float, default=1).fillna(1).astype(int)
        z = _get_first_existing(df, ['precursor_charge', 'charge', 'z'], cast=float, default=2).fillna(2).astype(int)

        rt = _get_first_existing(df, ['rt', 'retention_time', 'retention', 'retention_time_sec'], cast=float, default=np.nan)
        im = _get_first_existing(df, ['ion_mobility', 'mobility', 'ccs', 'k0'], cast=float, default=np.nan)

        pep_seq = df['peptide'].astype(str)
        proteins_raw = _get_first_existing(df, ['proteins', 'protein', 'protein_id'])
        proteins_raw = proteins_raw.astype(str) if proteins_raw is not None else pd.Series([''] * len(df))
        protein_ids, gene_ids, num_prot = self._split_accessions_and_entries(proteins_raw)

        # decoy detection from label 
        # Sage: label == -1 (decoy), +1 (target)
        label_series = pd.to_numeric(df['label'], errors='coerce')
        decoy = (label_series == -1)

        # spectrum-level q-value, peptide-level q-value and protein-level q-value
        pep = pd.to_numeric(df['posterior_error'], errors='coerce') if 'posterior_error' in df.columns else pd.Series([np.nan]*len(df))
        spectrum_q = pd.to_numeric(df['spectrum_q'], errors='coerce') if 'spectrum_q' in df.columns else pd.Series([np.nan]*len(df))
        peptide_q  = pd.to_numeric(df['peptide_q'], errors='coerce') if 'peptide_q' in df.columns else pd.Series([np.nan]*len(df))
        protein_q  = pd.to_numeric(df['protein_q'], errors='coerce') if 'protein_q' in df.columns else pd.Series([np.nan]*len(df))

        # precursor m/z if present (for joining with fragments)
        prec_mz = _get_first_existing(df, ['calcmass', 'expmass', 'precursor_mz', 'mz'], cast=float, default=np.nan)

        # modified peptide
        modpep = pep_seq.apply(self._annotate_unimod)

        # group id (same style as convert paths)
        group_id = run_id + "_" + scan_id.astype(str) + np.where(hit_rank > 1, "_rank" + hit_rank.astype(str), "")

        out = pd.DataFrame({
            'run_id': run_id,
            'scan_id': scan_id,
            'hit_rank': hit_rank,
            'massdiff': 0.0,
            'precursor_charge': z,
            'retention_time': rt,
            'ion_mobility': im,
            'peptide_sequence': pep_seq.str.replace(r'\[[-+0-9.]+\]', '', regex=True),
            'protein_id': protein_ids.fillna(''),
            'gene_id': gene_ids.fillna(''),
            'num_tot_proteins': num_prot.fillna(0).astype(int),
            'decoy': decoy.astype(bool),
            'modified_peptide': modpep,
            'group_id': group_id,
            'precursor_mz': prec_mz,
            'pep': pep,
            'q_value': spectrum_q,     
            'peptide_q': peptide_q,    
            'protein_q': protein_q     
        })
        return out


class SageFragmentParser:
    """
    Parse matched_fragments.sage.tsv → EasyPQP 'peaks' table:
      columns: scan_id, modified_peptide, precursor_charge, precursor_mz, fragment, product_mz, intensity
    """
    def __init__(self, frags_tsv: str, mz_precision_digits: int = 6):
        self.frags_tsv = frags_tsv
        self.mz_precision_digits = mz_precision_digits

    @staticmethod
    def _ann(ftype: str, ord_: int, z: int) -> str:
        return f"{ftype}{ord_}^{z}"

    def parse(self, psms_with_psmid: pd.DataFrame) -> pd.DataFrame:
        fr = pd.read_csv(self.frags_tsv, sep='\t', dtype=str).fillna('')
        for c in ['psm_id', 'fragment_ordinals', 'fragment_charge', 'fragment_mz_calculated',
                  'fragment_mz_experimental', 'fragment_intensity']:
            if c in fr.columns:
                fr[c] = pd.to_numeric(fr[c], errors='coerce')
        if 'psm_id' not in fr.columns:
            raise ValueError("matched_fragments.sage.tsv must contain a 'psm_id' column.")
        
        fr['psm_id'] = fr['psm_id'].astype(str).str.strip()

        fr['fragment'] = fr.apply(
            lambda r: self._ann(str(r['fragment_type']), int(r['fragment_ordinals']), int(r['fragment_charge'])),
            axis=1
        )
        fr['product_mz'] = fr['fragment_mz_calculated']

        # join to PSMs
        join_cols = ['psm_id', 'scan_id', 'modified_peptide', 'precursor_mz', 'precursor_charge', 'run_id']
        j = fr.merge(psms_with_psmid[join_cols], on='psm_id', how='inner')

        peaks = j[['run_id', 'scan_id', 'modified_peptide', 'precursor_charge', 'precursor_mz',
                   'fragment', 'product_mz', 'fragment_intensity']].copy()
        peaks.rename(columns={'fragment_intensity': 'intensity'}, inplace=True)

        # per-PSM normalization to 10,000 (matches convert paths)
        peaks['intensity'] = peaks['intensity'].fillna(0.0)
        grp = peaks.groupby(['run_id', 'scan_id', 'modified_peptide', 'precursor_charge'], dropna=False)['intensity']
        denom = grp.transform(lambda x: np.nanmax(x.values) if len(x) else np.nan)
        peaks['intensity'] = (peaks['intensity'] / denom) * 10000.0
        peaks['intensity'] = peaks['intensity'].fillna(0.0)

        # round and de-duplicate (keep most intense per exact fragment/product_mz)
        peaks['product_mz'] = peaks['product_mz'].round(self.mz_precision_digits)
        peaks['precursor_mz'] = peaks['precursor_mz'].round(self.mz_precision_digits)
        peaks['intensity'] = peaks['intensity'].round(self.mz_precision_digits)

        peaks = (peaks
                 .groupby(['run_id', 'scan_id', 'modified_peptide', 'precursor_charge', 'precursor_mz',
                           'fragment', 'product_mz'], as_index=False)['intensity']
                 .max())
        return peaks


def convert_sage(results_tsv: str,
                 fragments_tsv: str,
                 unimod_xml: Optional[str],
                 max_delta_unimod: float = 0.02) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    High-level conversion: Sage TSVs -> (psms_df, peaks_df)
    """
    # Read raw to extract psm_id for joining
    timestamped_echo("Info: Reading Sage PSMs")
    raw_res = pd.read_csv(results_tsv, sep='\t', dtype=str)
    if 'psm_id' not in raw_res.columns:
        raise ValueError("results.sage.tsv must contain a 'psm_id' for joining with matched fragments.")
    
    raw_res['psm_id'] = raw_res['psm_id'].astype(str).str.strip()

    psms = SagePSMParser(results_tsv, unimod_xml, max_delta_unimod).parse()
    psms = raw_res[['psm_id']].join(psms)
    
    if psms.empty:
        raise ValueError("No PSMs were parsed from the provided results.sage.tsv file.")

    timestamped_echo("Info: Reading Sage matched fragment peaks")
    peaks = SageFragmentParser(fragments_tsv).parse(psms)
    
    if peaks.empty:
        raise ValueError("No fragment peaks were parsed from the provided matched_fragments.sage.tsv file.")

    # Trim to minimal schema expected by library.generate
    keep = [
        'run_id', 'scan_id', 'hit_rank', 'massdiff', 'precursor_charge', 'retention_time',
        'ion_mobility', 'peptide_sequence', 'protein_id', 'gene_id', 'num_tot_proteins',
        'decoy', 'modified_peptide', 'group_id',
        'pep', 'q_value','peptide_q','protein_q'
    ]
    psms_export = psms[keep].copy()
    
    runs = sorted(psms_export['run_id'].dropna().unique().tolist())
    new_infiles = []
    for run in runs:
        psms_r  = psms_export.loc[psms_export['run_id'] == run]
        peaks_r = peaks.loc[peaks['run_id'] == run] if 'run_id' in peaks.columns else peaks

        if psms_r.empty or peaks_r.empty:
            timestamped_echo(f"Info: Skipping run {run}: psms={len(psms_r)}, peaks={len(peaks_r)}")
            continue

        psmpkl  = f"{run}.psmpkl"
        peakpkl = f"{run}.peakpkl"
        psms_r.to_pickle(psmpkl)
        peaks_r.to_pickle(peakpkl)
        timestamped_echo(f"Info: Wrote {psmpkl} and {peakpkl}")
        new_infiles.extend([psmpkl, peakpkl])

    if len(new_infiles) == 0:
        raise click.ClickException("No non-empty runs detected after Sage conversion.")

    return new_infiles
