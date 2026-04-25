#!/usr/bin/env python3
"""
Simulate CD deconvolution to test whether protein domains fold independently.

Based on the algebraic test in Lakshmanan, Hull, Berry, Burg, Bothner,
McKenna & Agbandje-McKenna, Viruses 2022, 14(9):1922.
https://doi.org/10.3390/v14091922

Their Eq. 2 recovers one domain's mean residue ellipticity (MRE) from
the full-length measurement and the other domain(s):

    MRE_i = (MRE_total × N_total  −  Σ_{j≠i} MRE_j × N_j) / N_i

If domains fold independently their MRE is additive (weighted by residue
count) and the recovered spectrum matches the directly measured one.

Usage:
    python simulate_cd_deconv.py                        # VP1u demo
    python simulate_cd_deconv.py --pdb 4AKE.pdb \\
        --chain A --domains "CORE:1-29,60-121;LID:122-159;NMP:30-59"
    python simulate_cd_deconv.py --fetch 1CLL \\
        --chain A --domains "N-lobe:4-75;C-lobe:82-148"
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

__version__ = "0.1.0"

# -- Style -----------------------------------------------------------------

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9,
    "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 250, "savefig.bbox": "tight", "savefig.facecolor": "white",
    "legend.framealpha": 0.9, "legend.edgecolor": "#cccccc",
})

COLORS = ["#2980b9", "#27ae60", "#8e44ad", "#e67e22",
          "#16a085", "#c0392b", "#2c3e50", "#d4ac0d"]
C_TOTAL, C_DECONV = "#2c3e50", "#e74c3c"
C_COUP_TOTAL, C_COUP_DECONV = "#c0392b", "#d35400"

WAV = np.arange(190.0, 261.0, 1.0)

# -- Basis spectra (Greenfield & Fasman parametric shapes) -----------------

def _g(x, mu, sig, amp):
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)

def basis_helix(w=WAV):
    """α-helix CD: positive ~191 nm, minima ~208 and ~222 nm."""
    return _g(w, 191, 5.5, 75) + _g(w, 208, 5.5, -33) + _g(w, 222, 7.5, -33)

def basis_sheet(w=WAV):
    """β-sheet CD: positive ~198 nm, minimum ~216 nm."""
    return _g(w, 198, 6, 20) + _g(w, 216, 8, -17)

def basis_disorder(w=WAV):
    """Random coil CD: strong negative ~198 nm."""
    return _g(w, 198, 7, -38) + _g(w, 218, 14, 4)

def mix_spectrum(f_h, f_e, f_c, w=WAV):
    """Weighted sum of basis spectra from secondary-structure fractions."""
    return f_h * basis_helix(w) + f_e * basis_sheet(w) + f_c * basis_disorder(w)

# -- Experimental error model ----------------------------------------------

def add_error(spectrum, rng, noise_sigma=0.35):
    """Heteroscedastic Gaussian noise + baseline drift + concentration error."""
    n = len(spectrum)
    wav_factor = 1.0 + 2.5 * np.exp(-((WAV[:n] - 190) / 4.5) ** 2)
    noise = rng.normal(0, noise_sigma, n) * wav_factor
    x = np.linspace(-1, 1, n)
    baseline = 0.25 * (rng.normal(0, 0.4) + rng.normal(0, 0.6) * x
                       + rng.normal(0, 0.25) * x**2)
    scale = 1.0 + rng.normal(0, 0.025)
    return spectrum * scale + noise + baseline

# -- Eq. 2 deconvolution --------------------------------------------------

def deconvolve(mre_total, n_total, others, n_target):
    """Recover one domain's MRE from full-length and all other domains (Eq. 2)."""
    return (mre_total * n_total - sum(m * n for m, n in others)) / n_target

# -- PDB parsing -----------------------------------------------------------

def parse_pdb_ss(path: str, chain: str = "A") -> Dict[int, str]:
    """Parse HELIX/SHEET records → {resnum: 'H'|'E'|'C'} for one chain."""
    resnums, helices, sheets = set(), [], []
    with open(path) as fh:
        for line in fh:
            if len(line) < 38:
                continue
            rec = line[:6].strip()
            if rec == "ATOM" and line[21] == chain:
                resnums.add(int(line[22:26]))
            elif rec == "HELIX" and line[19] == chain:
                helices.append((int(line[21:25]), int(line[33:37])))
            elif rec == "SHEET" and line[21] == chain:
                sheets.append((int(line[22:26]), int(line[33:37])))
    ss = {r: "C" for r in sorted(resnums)}
    for spans, tag in [(helices, "H"), (sheets, "E")]:
        for s, e in spans:
            for r in range(s, e + 1):
                if r in ss:
                    ss[r] = tag
    return ss

def fetch_pdb(pdb_id: str, dest: str) -> str:
    """Download PDB from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"Downloading {url} …")
    urllib.request.urlretrieve(url, dest)
    return dest

# -- Domain definitions ----------------------------------------------------

@dataclass
class DomainDef:
    name: str
    ranges: List[Tuple[int, int]]

    @property
    def all_residues(self) -> List[int]:
        return [r for s, e in self.ranges for r in range(s, e + 1)]

def parse_domain_string(s: str) -> List[DomainDef]:
    """Parse 'D1:1-90;D2:91-227' or 'CORE:1-29,60-121;LID:122-159'."""
    domains = []
    for part in s.split(";"):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"Bad domain spec '{part}' — expected NAME:start-end")
        name, ranges_str = part.split(":", 1)
        ranges = []
        for rng in ranges_str.split(","):
            toks = rng.strip().split("-")
            if len(toks) != 2:
                raise ValueError(f"Bad range '{rng}' in domain '{name}'")
            ranges.append((int(toks[0]), int(toks[1])))
        domains.append(DomainDef(name.strip(), ranges))
    if not domains:
        raise ValueError("No domains parsed — use format 'D1:1-90;D2:91-227'")
    return domains

def _ss_fracs(ss_map, residues):
    vals = [ss_map[r] for r in residues if r in ss_map]
    n = len(vals) or 1
    fh = sum(v == "H" for v in vals) / n
    fe = sum(v == "E" for v in vals) / n
    return fh, fe, 1.0 - fh - fe

# -- Simulation engine -----------------------------------------------------

@dataclass
class Domain:
    name: str
    n_res: int
    f_h: float
    f_e: float
    f_c: float
    mre_clean: np.ndarray = field(repr=False)
    mre_noisy: np.ndarray = field(repr=False)

@dataclass
class Scenario:
    label: str
    domains: List[Domain]
    n_total: int
    mre_total_noisy: np.ndarray
    deconv_noisy: List[np.ndarray]
    residuals: List[np.ndarray]
    perturbation: Optional[np.ndarray] = None

DomainSpec = Tuple[str, int, float, float, float]

def simulate(
    specs: List[DomainSpec],
    scenario: str = "independent",
    coupling_amp: float = 1.0,
    noise_sigma: float = 0.35,
    seed: int = 42,
) -> Scenario:
    rng = np.random.default_rng(seed)
    domains = []
    for name, n, fh, fe, fc in specs:
        clean = mix_spectrum(fh, fe, fc)
        domains.append(Domain(name, n, fh, fe, fc, clean,
                              add_error(clean, rng, noise_sigma)))
    n_total = sum(d.n_res for d in domains)
    if n_total == 0:
        raise ValueError("Total residue count is 0")

    mre_clean = sum(d.mre_clean * d.n_res for d in domains) / n_total
    pert = None
    if scenario == "coupled":
        pert = coupling_amp * (_g(WAV, 215, 8, -8) + _g(WAV, 232, 6, 5))
        mre_clean = mre_clean + pert
    mre_noisy = add_error(mre_clean, rng, noise_sigma)

    deconv, resid = [], []
    for i, d in enumerate(domains):
        others = [(domains[j].mre_noisy, domains[j].n_res)
                  for j in range(len(domains)) if j != i]
        rec = deconvolve(mre_noisy, n_total, others, d.n_res)
        deconv.append(rec)
        resid.append(rec - d.mre_noisy)

    return Scenario(scenario, domains, n_total, mre_noisy,
                    deconv, resid, pert)

# -- Figures ---------------------------------------------------------------

def _col(i):
    return COLORS[i % len(COLORS)]

def _style_ax(ax, xlabel="Wavelength (nm)", ylabel=None, title=None, title_color="k"):
    ax.axhline(0, lw=0.4, color="k")
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontweight="bold", fontsize=9, color=title_color)
    ax.legend(fontsize=6.5, frameon=True)


def _scenario_figure(sc, title_prefix, outpath):
    nd = len(sc.domains)
    fig, axes = plt.subplots(1, nd + 2, figsize=(4.2 * (nd + 2), 3.8))

    ax = axes[0]
    ct = C_TOTAL if sc.label == "independent" else C_COUP_TOTAL
    ax.plot(WAV, sc.mre_total_noisy, lw=1.8, color=ct, label="Full-length")
    for i, d in enumerate(sc.domains):
        ax.plot(WAV, d.mre_noisy, lw=1.4, color=_col(i), label=d.name)
    _style_ax(ax, ylabel="MRE (10³ deg·cm²·dmol⁻¹)", title="Measured spectra")

    ok = sc.label == "independent"
    cd = C_DECONV if ok else C_COUP_DECONV
    for i, d in enumerate(sc.domains):
        ax = axes[i + 1]
        ax.plot(WAV, d.mre_noisy, lw=2, color=_col(i), label=f"Measured {d.name}")
        ax.plot(WAV, sc.deconv_noisy[i], lw=1.5, ls="--", color=cd,
                label="Deconvolved (Eq. 2)")
        suffix = "" if ok else " — FAILS"
        _style_ax(ax, title=f"Recover {d.name}{suffix}",
                  title_color="k" if ok else "#c0392b")

    ax = axes[-1]
    for i, d in enumerate(sc.domains):
        ax.plot(WAV, sc.residuals[i], lw=1.3, color=_col(i), label=d.name)
    ax.fill_between(WAV, -1.5, 1.5, alpha=0.07, color="gray", label="Noise floor")
    _style_ax(ax, ylabel="Δ MRE", title="Residuals (deconv − measured)")

    fig.suptitle(title_prefix, fontsize=11, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def _comparison_figure(indep, coupled, outpath, protein_name):
    d0 = indep.domains[0].name
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))

    for ax, sc, cd, title, tc in [
        (axes[0], indep, C_DECONV, "Independent: close match", "#1a5276"),
        (axes[1], coupled, C_COUP_DECONV, "Coupled: systematic mismatch", "#c0392b"),
    ]:
        ax.plot(WAV, sc.domains[0].mre_noisy, lw=2.2, color=COLORS[0],
                label=f"Measured {d0}")
        ax.plot(WAV, sc.deconv_noisy[0], lw=1.8, ls="--", color=cd,
                label="Deconvolved (Eq. 2)")
        _style_ax(ax, ylabel="MRE", title=title, title_color=tc)

    ax = axes[2]
    ax.plot(WAV, indep.residuals[0], lw=1.8, color=C_TOTAL, label="Independent")
    ax.plot(WAV, coupled.residuals[0], lw=1.8, color=C_COUP_TOTAL, label="Coupled")
    ax.fill_between(WAV, -1.5, 1.5, alpha=0.08, color="gray", label="Noise floor")
    _style_ax(ax, ylabel="Δ MRE", title="Residuals reveal coupling")

    fig.suptitle(f"Eq. 2 deconvolution test — {protein_name}\n"
                 f"(simulated data, {d0} recovery shown)",
                 fontsize=11, fontweight="bold", y=1.04)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def _basis_figure(outpath):
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    for fn, label, c in [(basis_helix, "α-helix", "#e74c3c"),
                          (basis_sheet, "β-sheet", "#3498db"),
                          (basis_disorder, "Disorder / coil", "#95a5a6")]:
        ax.plot(WAV, fn(), lw=2, color=c, label=label)
    _style_ax(ax, ylabel="Δε (10³ deg·cm²·dmol⁻¹)",
              title="Basis CD spectra (Greenfield & Fasman shapes)")
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def generate_figures(specs, outdir, protein_name="protein",
                     coupling_amp=1.0, noise_sigma=0.35, seed=42):
    os.makedirs(outdir, exist_ok=True)
    _basis_figure(os.path.join(outdir, "basis_spectra.png"))

    indep = simulate(specs, "independent", noise_sigma=noise_sigma, seed=seed)
    _scenario_figure(indep, f"Scenario A — Independent folding ({protein_name})",
                     os.path.join(outdir, "independent_deconv.png"))

    coupled = simulate(specs, "coupled", coupling_amp=coupling_amp,
                       noise_sigma=noise_sigma, seed=seed + 1)
    _scenario_figure(coupled, f"Scenario B — Coupled folding ({protein_name})",
                     os.path.join(outdir, "coupled_deconv.png"))

    _comparison_figure(indep, coupled, os.path.join(outdir, "comparison.png"),
                       protein_name)

    for f in ["basis_spectra", "independent_deconv", "coupled_deconv", "comparison"]:
        print(f"  {f}.png")
    print(f"\nAll figures → {outdir}/")

# -- VP1u demo (Table 1 from Lakshmanan et al. 2022) -----------------------

VP1U_SPECS: List[DomainSpec] = [
    ("RBD (D1)",  90,  0.40, 0.00, 0.60),
    ("PLA2 (D2)", 137, 0.10, 0.25, 0.65),
]

def run_demo(outdir):
    print("─── B19V VP1u demo (Lakshmanan et al. 2022, Table 1) ───")
    n_tot = sum(n for _, n, *_ in VP1U_SPECS)
    print(f"  VP1u full-length: {n_tot} aa")
    for name, n, fh, fe, fc in VP1U_SPECS:
        print(f"  {name}: {n} aa  (H={fh:.0%} E={fe:.0%} C={fc:.0%})")
    print()
    generate_figures(VP1U_SPECS, outdir, protein_name="B19V VP1u")

# -- PDB pipeline ----------------------------------------------------------

def run_pdb(pdb_path, chain, domain_defs, outdir, coupling_amp=1.0, noise_sigma=0.35):
    if not os.path.exists(pdb_path):
        print(f"ERROR: {pdb_path} not found"); sys.exit(1)
    ss_map = parse_pdb_ss(pdb_path, chain)
    if not ss_map:
        print(f"ERROR: no ATOM records for chain {chain} in {pdb_path}"); sys.exit(1)

    n_all = len(ss_map)
    n_h = sum(v == "H" for v in ss_map.values())
    n_e = sum(v == "E" for v in ss_map.values())
    print(f"─── PDB: {pdb_path}  chain {chain} ───")
    print(f"  {n_all} residues  H={n_h} ({n_h/n_all:.0%})  "
          f"E={n_e} ({n_e/n_all:.0%})  C={n_all-n_h-n_e} ({(n_all-n_h-n_e)/n_all:.0%})")

    specs = []
    for dd in domain_defs:
        residues = [r for r in dd.all_residues if r in ss_map]
        n = len(residues)
        if n == 0:
            print(f"  WARNING: {dd.name} has 0 residues in chain {chain}, skipping")
            continue
        fh, fe, fc = _ss_fracs(ss_map, residues)
        specs.append((dd.name, n, fh, fe, fc))
        print(f"  {dd.name}: {n} aa  (H={fh:.0%} E={fe:.0%} C={fc:.0%})")

    if not specs:
        print("ERROR: no valid domains"); sys.exit(1)

    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0].upper()
    print()
    generate_figures(specs, outdir, protein_name=f"{pdb_name} chain {chain}",
                     coupling_amp=coupling_amp, noise_sigma=noise_sigma)

# -- CLI -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="CD deconvolution simulation for domain independence testing.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  %(prog)s                                          # VP1u demo
  %(prog)s --fetch 1CLL --chain A \\
           --domains "N-lobe:4-75;C-lobe:82-148"    # auto-download from RCSB
  %(prog)s --pdb my.pdb --chain A \\
           --domains "CORE:1-29,60-121;LID:122-159" # local PDB, discontinuous domains
""")
    ap.add_argument("--pdb", help="Path to a local PDB file")
    ap.add_argument("--fetch", metavar="PDB_ID",
                    help="4-letter PDB ID to download from RCSB (cached locally)")
    ap.add_argument("--chain", default="A", help="Chain ID (default: A)")
    ap.add_argument("--domains",
                    help="Domain definitions, e.g. 'D1:1-90;D2:91-227'")
    ap.add_argument("--outdir", default=None,
                    help="Output directory (default: figures/ or figures_pdb/)")
    ap.add_argument("--coupling", type=float, default=1.0,
                    help="Coupling perturbation amplitude (default: 1.0)")
    ap.add_argument("--noise", type=float, default=0.35,
                    help="Noise sigma on MRE scale (default: 0.35)")
    ap.add_argument("--seed", type=int, default=42, help="RNG seed (default: 42)")
    args = ap.parse_args()

    base = os.path.dirname(os.path.abspath(__file__))

    if args.pdb or args.fetch:
        if not args.domains:
            ap.error("--domains is required when using --pdb or --fetch")
        domain_defs = parse_domain_string(args.domains)
        if args.fetch:
            pdb_path = os.path.join(base, f"{args.fetch.upper()}.pdb")
            if not os.path.exists(pdb_path):
                fetch_pdb(args.fetch, pdb_path)
            else:
                print(f"Using cached {pdb_path}")
        else:
            pdb_path = args.pdb
        outdir = args.outdir or os.path.join(base, "figures_pdb")
        run_pdb(pdb_path, args.chain, domain_defs, outdir,
                coupling_amp=args.coupling, noise_sigma=args.noise)
    else:
        if args.domains:
            ap.error("--domains requires --pdb or --fetch")
        outdir = args.outdir or os.path.join(base, "figures")
        run_demo(outdir)

if __name__ == "__main__":
    main()
