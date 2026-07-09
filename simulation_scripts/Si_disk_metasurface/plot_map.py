"""
Plot the angle-resolved reflectance energy-momentum maps for the bare a-Si disk
BIC metasurface. Produces, for each of the four map_*.npz files, a two-panel PNG
(left: full reflectance range; right: contrast-enhanced to reveal the weak,
absorption-broadened BIC band) plus a 2x2 overview figure.

Usage: python plot_map.py [outbase_dir]
If no dir is given, uses the default Si_disk_20260708_bare-reflectance folder.
"""
import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT = os.path.join(REPO, "awesome_photonics_sim_output",
                       "Si_disk_20260708_bare-reflectance")

TITLES = {
    "1a_ypol_ky-sweep_kx0": "1a  y-pol, ky-sweep (kx=0)  [TM-like]",
    "1b_ypol_kx-sweep_ky0": "1b  y-pol, kx-sweep (ky=0)  [TE-like, BIC]",
    "2a_xpol_ky-sweep_kx0": "2a  x-pol, ky-sweep (kx=0)  [TE-like, BIC]",
    "2b_xpol_kx-sweep_ky0": "2b  x-pol, kx-sweep (ky=0)  [TM-like]",
}
ORDER = ["1a_ypol_ky-sweep_kx0", "1b_ypol_kx-sweep_ky0",
         "2a_xpol_ky-sweep_kx0", "2b_xpol_kx-sweep_ky0"]

BIC_VMAX = 0.06  # contrast-enhanced ceiling to reveal the damped BIC band


def _load(outbase, label):
    f = os.path.join(outbase, label, f"map_{label}.npz")
    if not os.path.exists(f):
        return None
    d = np.load(f, allow_pickle=True)
    return d["wavelength_nm"], d["angle_deg"], d["R"], f


def plot_single(outbase, label):
    got = _load(outbase, label)
    if got is None:
        print("missing", label)
        return
    wl, ang, R, f = got
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, (vmax, tag) in zip(axes, [(max(0.15, float(np.percentile(R, 99.5))), "full range"),
                                      (BIC_VMAX, f"enhanced (vmax={BIC_VMAX})")]):
        im = ax.pcolormesh(ang, wl, R, shading="auto", cmap="gray", vmin=0, vmax=vmax)
        ax.set_xlabel("Incidence angle (deg)")
        ax.set_ylabel("Wavelength (nm)")
        ax.set_title(tag, fontsize=10)
        fig.colorbar(im, ax=ax, label="Reflectance")
    fig.suptitle(TITLES.get(label, label), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(outbase, label, f"map_{label}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def _geom_label(outbase):
    for label in ORDER:
        f = os.path.join(outbase, label, f"map_{label}.npz")
        if os.path.exists(f):
            d = np.load(f, allow_pickle=True)
            if "period_nm" in d:
                return f"P={float(d['period_nm']):.0f}, R={float(d['radius_nm']):.0f}, H={float(d['height_nm']):.0f} nm"
    return ""


def plot_overview(outbase):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, label in zip(axes.ravel(), ORDER):
        got = _load(outbase, label)
        if got is None:
            ax.set_visible(False)
            continue
        wl, ang, R, f = got
        im = ax.pcolormesh(ang, wl, R, shading="auto", cmap="gray", vmin=0, vmax=BIC_VMAX)
        ax.set_title(TITLES.get(label, label), fontsize=9)
        ax.set_xlabel("Angle (deg)")
        ax.set_ylabel("Wavelength (nm)")
        fig.colorbar(im, ax=ax, label="R")
    fig.suptitle(f"Bare a-Si disk BIC metasurface ({_geom_label(outbase)}) "
                 f"- total reflectance, enhanced (vmax={BIC_VMAX})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(outbase, "overview_4maps_enhanced.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    outbase = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    for label in ORDER:
        plot_single(outbase, label)
    plot_overview(outbase)
