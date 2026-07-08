"""Plot a single energy-momentum map (wavelength vs angle) from map_<label>.npz."""
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

label = sys.argv[1] if len(sys.argv) > 1 else "kx"
fname = f"map_{label}.npz"
d = np.load(fname)
wl = d["wavelength_nm"]
ang = d["angle_deg"]
R = d["R"]

axis_name = {"kx": "angle in k$_x$ plane (⊥ wires, deg)", "ky": "angle in k$_y$ plane (∥ wires, deg)"}[label]
title = {
    "kx": "TM reflectance, k$_x$-sweep (k$_y$=0)",
    "ky": "TM reflectance, k$_y$-sweep (k$_x$=0)",
}[label]

fig, ax = plt.subplots(figsize=(7, 5.5))
pc = ax.pcolormesh(ang, wl, R, shading="auto", cmap="inferno", vmin=0, vmax=1)
ax.set_xlabel(axis_name)
ax.set_ylabel("wavelength (nm)")
ax.set_title(title)
cb = fig.colorbar(pc, ax=ax)
cb.set_label("Reflectance (specular)")
fig.tight_layout()
outname = f"map_{label}.png"
fig.savefig(outname, dpi=180)
print(f"saved {outname}  R range [{np.nanmin(R):.3f},{np.nanmax(R):.3f}]  shape={R.shape}  nan_count={np.isnan(R).sum()}")
