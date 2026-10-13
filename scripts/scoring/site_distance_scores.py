#!/usr/bin/env python3
"""
Site-specific ligand distances for co-folding models (AF3 / Boltz2 / RoseTTAFold-AllAtom).

Generalises the TGFBR2 SS1 scorer (TGFb_full/scripts/fold_ss1.py): the site residues, the
residues to report individually, the contact cutoff and the reference are all flags.

For each compound model:
  <site>_min_dist             closest ligand heavy atom to any site residue heavy atom (A)
  <site>_residues_contacted   site residues with a ligand heavy atom < --contact-cutoff
  <site>_n_residues           site residues found in the model (fewer = numbering/mapping gap)
  <site>_contacted            comma-separated contacted residues, e.g. "D32,K50"
  <site>_centroid_dist        ligand center-of-mass to site center-of-mass (A); distinguishes a
                              ligand actually sitting in the pocket from one atom grazing its edge
  res<N>_dist                 closest ligand heavy atom to residue N (one per --key-residues)
  <site>_rmsd                 site heavy-atom RMSD after a global CA fit (only with --reference)

Residue numbers are in the --reference numbering when given (mapped onto each model by
sequence alignment, so models with different numbering still work); otherwise in the
model's own numbering. Distances are in the model's own frame, so no superposition is needed.
Ligand = every non-polymer, non-water residue (AF3/Boltz2 chain Z, RFAA LG1, ...).

Usage:
  python site_distance_scores.py --method af3 --results-folder <...>/fine_screening/AF3/output \
      --site 30-32,50-55 --site-name ss1 --key-residues 32 \
      --reference Input/protein_file/TGFBR2/TGFBR2.pdb --output <...>/AF3/site_ss1_scores.csv
"""
import argparse
import glob
import os
import sys
from multiprocessing import Pool

import gemmi
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

# Per-method layout: glob for compound dirs, glob for the top-ranked model inside one
LAYOUT = {
    "af3": ("*", "*_model.cif"),
    "boltz2": ("job_*", "*_model_0.cif"),
    "rfaa": ("compound_*", "structure.pdb"),
}
WATER = {"HOH", "WAT", "DOD"}


def parse_residues(spec):
    """'30-32,50-55,D97' -> [30, 31, 32, 50, ..., 55, 97] (letters ignored, order kept, no dups)."""
    out = []
    for part in filter(None, (p.strip() for p in spec.split(","))):
        part = part.lstrip("ACDEFGHIKLMNPQRSTVWY")  # allow one-letter prefixes like D32
        lo, _, hi = part.partition("-")
        try:
            nums = range(int(lo), int(hi or lo) + 1)
        except ValueError:
            sys.exit(f"ERROR: bad residue spec '{part}' in '{spec}' (use e.g. 30-32,50-55)")
        out += [n for n in nums if n not in out]
    return out


def protein_residues(chain):
    """Residues with a CA atom, plus their one-letter sequence."""
    res = [r for r in chain if r.find_atom("CA", "*") is not None]
    seq = "".join((gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                   if gemmi.find_tabulated_residue(r.name) else "X") or "X" for r in res)
    return seq, res


def pick_chain(model, name=None):
    if name:
        for ch in model:
            if ch.name == name:
                return ch
        raise ValueError(f"chain {name} not found (have {[c.name for c in model]})")
    return max(model, key=lambda ch: len(protein_residues(ch)[1]))  # receptor = most residues


def is_polymer(res):
    info = gemmi.find_tabulated_residue(res.name)
    return info is not None and (info.is_amino_acid() or info.is_nucleic_acid())


def ligand_xyz(model):
    """Heavy atoms of the ligand: the largest non-polymer, non-water residue (by atom count).

    Picking only the largest residue -- rather than every hetero residue -- drops disconnected
    counterions (Cl-, Br-, Na+, K+, ...) that ship in the same SMILES (e.g. "R.Cl"); a lone ion
    can land anywhere in the model and would otherwise skew centroid_dist/min_dist/contacts
    toward wherever it happens to sit, not toward the actual drug-like ligand.
    """
    residues = [r for ch in model for r in ch if r.name not in WATER and not is_polymer(r)]
    if not residues:
        return np.empty((0, 3))
    lig = max(residues, key=lambda r: sum(1 for a in r if a.element.name != "H"))
    return xyz(lig)


def xyz(atoms):
    return np.array([[a.pos.x, a.pos.y, a.pos.z] for a in atoms if a.element.name != "H"]).reshape(-1, 3)


def superpose(a, b):
    """Kabsch fit of a onto b; returns a function applying it."""
    ma, mb = a.mean(0), b.mean(0)
    u, _, vt = np.linalg.svd((a - ma).T @ (b - mb))
    rot = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
    return lambda x: (x - ma) @ rot + mb


def min_dist(lig, pts):
    if len(lig) == 0 or len(pts) == 0:
        return np.nan
    return float(np.sqrt(((lig[:, None] - pts[None]) ** 2).sum(-1)).min())


def load_reference(path, chain):
    st = gemmi.read_structure(path)
    st.remove_hydrogens()
    return protein_residues(pick_chain(st[0], chain))


def map_to_reference(mseq, mres, ref):
    """Align model receptor to reference; return {ref_resnum: model_residue} and CA pairs."""
    ref_seq, ref_res = ref
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    aln = aligner.align(mseq, ref_seq)[0]
    pairs = [(i, j) for (a0, a1), (b0, b1) in zip(*aln.aligned)
             for i, j in zip(range(a0, a1), range(b0, b1))]
    return {ref_res[j].seqid.num: mres[i] for i, j in pairs}, pairs


def init_worker(opts):
    """Each worker loads its own reference: gemmi objects can't be pickled between processes."""
    global OPTS, REF
    OPTS = opts
    REF = load_reference(opts.reference, opts.reference_chain) if opts.reference else None


def score(cdir):
    opts, ref = OPTS, REF
    row = {"folder": os.path.basename(cdir)}
    smi = os.path.join(cdir, "smiles.txt")
    if os.path.isfile(smi):
        row["SMILES"] = open(smi).read().strip()
    s = opts.site_name
    try:
        st = gemmi.read_structure(sorted(glob.glob(os.path.join(cdir, opts.model_glob)))[0])
        st.remove_hydrogens()
        model = st[0]
        rec = pick_chain(model, opts.receptor_chain)
        mseq, mres = protein_residues(rec)
        if ref:
            by_num, pairs = map_to_reference(mseq, mres, ref)
        else:
            by_num, pairs = {r.seqid.num: r for r in mres}, None
        lig = ligand_xyz(model)
        if len(lig) == 0:
            raise ValueError("no ligand atoms found")

        site = {n: by_num[n] for n in opts.site if n in by_num}
        d = {n: min_dist(lig, xyz(r)) for n, r in site.items()}
        contacted = [f"{gemmi.find_tabulated_residue(site[n].name).one_letter_code.upper()}{n}"
                     for n in site if d[n] < opts.contact_cutoff]
        row[f"{s}_min_dist"] = round(min(d.values()), 2) if d else np.nan
        row[f"{s}_residues_contacted"] = len(contacted)
        row[f"{s}_n_residues"] = len(site)
        row[f"{s}_contacted"] = ",".join(contacted)
        if site:
            site_xyz = np.concatenate([xyz(r) for r in site.values()])
            row[f"{s}_centroid_dist"] = round(float(np.linalg.norm(lig.mean(0) - site_xyz.mean(0))), 2)
        for n in opts.key_residues:
            row[f"res{n}_dist"] = round(min_dist(lig, xyz(by_num[n])), 2) if n in by_num else np.nan

        if ref:  # site RMSD after a global CA fit onto the reference
            ref_res = ref[1]
            fit = superpose(xyz(mres[i]["CA"][0] for i, _ in pairs), xyz(ref_res[j]["CA"][0] for _, j in pairs))
            mp, rp = [], []
            for i, j in pairs:
                if ref_res[j].seqid.num in opts.site:
                    for a in ref_res[j]:
                        b = mres[i].find_atom(a.name, "*")
                        if b is not None and a.element.name != "H":
                            mp.append(b)
                            rp.append(a)
            if mp:
                row[f"{s}_rmsd"] = round(float(np.sqrt(((fit(xyz(mp)) - xyz(rp)) ** 2).sum(1).mean())), 3)
    except Exception as exc:
        row["error"] = repr(exc)[:200]
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--method", choices=LAYOUT, required=True)
    parser.add_argument("--results-folder", required=True,
                        help="Method output root containing one folder per compound (job_* / compound_*)")
    parser.add_argument("--site", required=True, help="Site residues, e.g. 30-32,50-55 (or D32,K50)")
    parser.add_argument("--site-name", default="site", help="Column prefix, e.g. ss1")
    parser.add_argument("--key-residues", default="", help="Residues reported individually as res<N>_dist")
    parser.add_argument("--contact-cutoff", type=float, default=4.0, help="Contact distance (A)")
    parser.add_argument("--reference", default=None,
                        help="Reference receptor PDB/CIF; residue numbers then refer to it (+ adds <site>_rmsd)")
    parser.add_argument("--reference-chain", default=None, help="Reference chain (default: longest)")
    parser.add_argument("--receptor-chain", default=None, help="Model receptor chain (default: longest)")
    parser.add_argument("--model-glob", default=None, help="Override model file glob inside each folder")
    parser.add_argument("--output", required=True, help="Output CSV")
    parser.add_argument("--nproc", type=int, default=int(os.environ.get("NPROC", 4)))
    parser.add_argument("--limit", type=int, default=None, help="Score only the first N folders (testing)")
    opts = parser.parse_args()

    opts.site = parse_residues(opts.site)
    opts.key_residues = parse_residues(opts.key_residues) if opts.key_residues else []
    dir_glob, default_model = LAYOUT[opts.method]
    opts.model_glob = opts.model_glob or default_model
    ref = load_reference(opts.reference, opts.reference_chain) if opts.reference else None
    if ref:
        missing = [n for n in opts.site + opts.key_residues if n not in {r.seqid.num for r in ref[1]}]
        if missing:
            print(f"WARNING: residues {missing} not in reference {opts.reference}")

    dirs = sorted(d for d in glob.glob(os.path.join(opts.results_folder, dir_glob))
                  if os.path.isdir(d) and glob.glob(os.path.join(d, opts.model_glob)))
    if opts.limit:
        dirs = dirs[: opts.limit]
    if not dirs:
        sys.exit(f"ERROR: no {opts.method} models ({dir_glob}/{opts.model_glob}) under {opts.results_folder}")
    print(f"{opts.method}: {len(dirs)} models; site {opts.site_name} = {opts.site}; "
          f"key {opts.key_residues}; cutoff {opts.contact_cutoff} A; reference {opts.reference or 'none'}")

    with Pool(opts.nproc, initializer=init_worker, initargs=(opts,)) as pool:
        rows = pool.map(score, dirs, chunksize=20)
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(os.path.abspath(opts.output)), exist_ok=True)
    df.to_csv(opts.output, index=False)

    s = opts.site_name
    n_err = int(df["error"].notna().sum()) if "error" in df else 0
    print(f"Saved: {opts.output}  ({len(df)} rows, {n_err} errors)")
    if f"{s}_residues_contacted" in df:
        print(f"contacting {s} (>=1 residue < {opts.contact_cutoff} A): {int((df[f'{s}_residues_contacted'] > 0).sum())}")
        num = [c for c in df if c.endswith(("_dist", "_rmsd", "_contacted")) and df[c].dtype != object]
        print(df[num].describe(percentiles=[.1, .5, .9]).round(2).to_string())


if __name__ == "__main__":
    main()
