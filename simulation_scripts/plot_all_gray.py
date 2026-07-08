"""
Generate a single-panel GRAYSCALE reflectance energy-momentum PNG for every
map_*.npz across all simulation-output projects (AlGaAs, SiN, Si-disk).

Each npz shares the schema: wavelength_nm (nwl,), angle_deg (nang,), R (nwl,nang).
Contrast is auto-scaled per map (vmin=1st pct, vmax=99th pct) so both the strong
(AlGaAs/SiN) and the weak (lossy Si-disk BIC) resonances are visible.

Output: <npz_stem>_gray.png next to each npz (existing colour PNGs untouched).

Run: python plot_all_gray.py
"""
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

GAMMA = 0.5  # sqrt stretch: brightens weak features (e.g. the lossy Si-disk BIC)

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUTROOT = os.path.join(REPO, "awesome_photonics_sim_output")


def plot_one(npz_path):
    d = np.load(npz_path, allow_pickle=True)
    if "R" not in d or "wavelength_nm" not in d or "angle_deg" not in d:
        return None
    wl = d["wavelength_nm"]
    ang = d["angle_deg"]
    R = np.asarray(d["R"], dtype=float)
    finite = R[np.isfinite(R)]
    vmin = float(np.percentile(finite, 1))
    vmax = float(np.percentile(finite, 99))
    if vmax <= vmin:
        vmax = vmin + 1e-6

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    norm = PowerNorm(gamma=GAMMA, vmin=max(vmin, 0.0), vmax=vmax)
    im = ax.pcolormesh(ang, wl, R, shading="auto", cmap="gray", norm=norm)
    ax.set_xlabel("Incidence angle (deg)")
    ax.set_ylabel("Wavelength (nm)")
    rel = os.path.relpath(npz_path, OUTROOT).replace("\\", "/")
    ax.set_title(rel, fontsize=8)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(f"Reflectance  (grayscale, gamma={GAMMA}, max {vmax:.3f})")
    fig.tight_layout()
    out = os.path.splitext(npz_path)[0] + "_gray.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


if __name__ == "__main__":
    npzs = sorted(glob.glob(os.path.join(OUTROOT, "**", "*.npz"), recursive=True))
    n = 0
    for f in npzs:
        out = plot_one(f)
        if out:
            n += 1
            print("wrote", os.path.relpath(out, OUTROOT))
    print(f"\n{n} grayscale maps written under {OUTROOT}")
