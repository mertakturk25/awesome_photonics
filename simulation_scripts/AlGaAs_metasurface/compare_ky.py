"""
Quantitative comparison of the simulated TM ky energy-momentum map
(map_ky.npz, now full 500-900nm x +-60deg) against the experimental
grayscale image "ALGAS_PHOTO_1 - Copy.jpg" (image cropped to exactly
+-60deg in x; wavelength axis spans 500-900nm top(900)-to-bottom(500)).
Both now cover the identical range -- no sub-cropping needed.
"""
import numpy as np
from PIL import Image
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

IMG_PATH = r"C:\Users\merta\Documents\Github\awesome_photonics\ALGAS_PHOTO_1 - Copy.jpg"
WL_TOP, WL_BOT = 900.0, 500.0  # experimental image wavelength axis (top->bottom)
ANGLE_MIN, ANGLE_MAX = -60.0, 60.0

img = np.array(Image.open(IMG_PATH).convert("L")).astype(float)
H, W = img.shape
angle_img = np.linspace(ANGLE_MIN, ANGLE_MAX, W)
wl_img = np.linspace(WL_TOP, WL_BOT, H)  # decreasing (row0=900nm ... rowH-1=500nm)

# RegularGridInterpolator needs increasing coordinates -> flip rows
wl_img_inc = wl_img[::-1]
img_inc = img[::-1, :]
interp_img = RegularGridInterpolator((wl_img_inc, angle_img), img_inc, bounds_error=False, fill_value=np.nan)

d = np.load("map_ky.npz")
wl_sim = d["wavelength_nm"]
angle_sim = d["angle_deg"]
R_sim = d["R"]

AA, WW = np.meshgrid(angle_sim, wl_sim)
img_on_sim_grid = interp_img((WW, AA))

# normalize both to [0,1] for shape comparison (camera intensity is not
# calibrated absolute reflectance, so only relative/normalized comparison
# is meaningful)
def norm(x):
    x = x - np.nanmin(x)
    return x / np.nanmax(x)

img_n = norm(img_on_sim_grid)
sim_n = norm(R_sim)

n_bad_sim = int(np.isnan(R_sim).sum())
if n_bad_sim:
    print(f"note: {n_bad_sim} NaN point(s) in simulated R (isolated RCWA solver failure); excluded from stats")
mask = ~np.isnan(img_n) & ~np.isnan(sim_n)
pearson_r = np.corrcoef(img_n[mask].ravel(), sim_n[mask].ravel())[0, 1]
rmse = np.sqrt(np.nanmean((img_n - sim_n) ** 2))

# ridge (peak-intensity wavelength per angle) comparison
ridge_img = wl_sim[np.nanargmax(img_n, axis=0)]
ridge_sim = wl_sim[np.nanargmax(sim_n, axis=0)]
ridge_diff = ridge_img - ridge_sim

print(f"Pearson correlation (normalized): {pearson_r:.4f}")
print(f"RMSE (normalized 0-1 scale): {rmse:.4f}")
print(f"Ridge (peak-lambda vs angle) mean abs diff: {np.mean(np.abs(ridge_diff)):.1f} nm, "
      f"median: {np.median(np.abs(ridge_diff)):.1f} nm, max: {np.max(np.abs(ridge_diff)):.1f} nm")

fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
im0 = axes[0].pcolormesh(angle_sim, wl_sim, img_n, shading="auto", cmap="gray", vmin=0, vmax=1)
axes[0].set_title("Experimental (cropped, normalized)")
axes[0].set_xlabel("angle (deg)")
axes[0].set_ylabel("wavelength (nm)")
fig.colorbar(im0, ax=axes[0])

im1 = axes[1].pcolormesh(angle_sim, wl_sim, sim_n, shading="auto", cmap="inferno", vmin=0, vmax=1)
axes[1].set_title("Simulated TM ky-sweep (normalized)")
axes[1].set_xlabel("angle (deg)")
fig.colorbar(im1, ax=axes[1])

im2 = axes[2].pcolormesh(angle_sim, wl_sim, img_n - sim_n, shading="auto", cmap="coolwarm", vmin=-1, vmax=1)
axes[2].set_title(f"Difference (exp-sim)\nr={pearson_r:.3f}, RMSE={rmse:.3f}")
axes[2].set_xlabel("angle (deg)")
fig.colorbar(im2, ax=axes[2])

fig.tight_layout()
fig.savefig("compare_ky.png", dpi=170)
print("saved compare_ky.png")

fig2, ax2 = plt.subplots(figsize=(6.5, 5))
ax2.plot(angle_sim, ridge_img, label="experimental peak-intensity wavelength")
ax2.plot(angle_sim, ridge_sim, label="simulated peak-R wavelength")
ax2.set_xlabel("angle (deg)")
ax2.set_ylabel("wavelength of brightest/highest-R pixel (nm)")
ax2.legend()
ax2.set_title("Ridge (dispersion branch) comparison, 650-800nm window")
fig2.tight_layout()
fig2.savefig("compare_ky_ridge.png", dpi=170)
print("saved compare_ky_ridge.png")
