"""
Compute the four angle-resolved TOTAL-reflectance energy-momentum maps for the
bare a-Si disk BIC metasurface (Marangi et al., arXiv:2511.08103, Fig. 2a --
no J-aggregate coating), lambda 500-800 nm x theta [-40, +40] deg:

  1a  y-pol, ky-sweep (kx=0)   phi=90, psi=0    (TM-like: BIC dark)
  1b  y-pol, kx-sweep (ky=0)   phi=0,  psi=90   (TE-like: BIC bright)
  2a  x-pol, ky-sweep (kx=0)   phi=90, psi=90   (TE-like: BIC bright)
  2b  x-pol, kx-sweep (ky=0)   phi=0,  psi=0    (TM-like: BIC dark)

For an ideal circular disk on a square lattice (C4v symmetry) a 90-deg rotation
maps x<->y and kx<->ky, so exactly:
      2a == 1b   and   2b == 1a
This was verified numerically to ~1e-14 (see README). We therefore compute only
the two distinct y-polarized maps (1a, 1b) and write 2a, 2b via this exact
symmetry relation. (A real, slightly asymmetric / C2 pillar would break the
degeneracy and make the cross-polarization maps independent.)

Each map is saved into its own subfolder under
  awesome_photonics_sim_output/Si_disk_20260708_bare-reflectance/<label>/

Run: python sweep.py [nh] [wl_step] [angle_step]
"""
import os
import sys
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import total_R, P, R, H

NH = int(sys.argv[1]) if len(sys.argv) > 1 else 180
WL_STEP = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
ANGLE_STEP = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
DISC = (256, 256)

WL_MIN, WL_MAX = 500.0, 800.0
ANG_MIN, ANG_MAX = -40.0, 40.0
WLS = np.arange(WL_MIN, WL_MAX + 1e-6, WL_STEP)
ANGLES = np.arange(ANG_MIN, ANG_MAX + 1e-6, ANGLE_STEP)

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUTBASE = os.path.join(REPO, "awesome_photonics_sim_output",
                       "Si_disk_20260708_bare-reflectance")

# (phi, psi) for the two distinct computed maps (both y-polarized)
COMPUTED = {
    "1a_ypol_ky-sweep_kx0": dict(phi=90.0, psi=0.0, note="TM-like (BIC dark)"),
    "1b_ypol_kx-sweep_ky0": dict(phi=0.0, psi=90.0, note="TE-like (BIC bright)"),
}
# symmetry-mirrored maps: label -> (source label, phi, psi it is equivalent to)
MIRRORED = {
    "2a_xpol_ky-sweep_kx0": dict(source="1b_ypol_kx-sweep_ky0", phi=90.0, psi=90.0),
    "2b_xpol_kx-sweep_ky0": dict(source="1a_ypol_ky-sweep_kx0", phi=0.0, psi=0.0),
}


def _worker(args):
    wl, theta, phi, psi = args
    return total_R(wl, theta, phi, psi, nh=NH, disc=DISC)


def run_one(label, phi, psi, note):
    outdir = os.path.join(OUTBASE, label)
    os.makedirs(outdir, exist_ok=True)
    jobs = [(wl, th, phi, psi) for wl in WLS for th in ANGLES]
    n = len(jobs)
    print(f"[{label}] {n} pts (nh={NH}, {len(WLS)}wl x {len(ANGLES)}ang), phi={phi} psi={psi} {note}", flush=True)
    t0 = time.time()
    results = np.zeros(n)
    done = 0
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        futs = {ex.submit(_worker, j): i for i, j in enumerate(jobs)}
        for fut in as_completed(futs):
            results[futs[fut]] = fut.result()
            done += 1
            if done % 1000 == 0:
                el = time.time() - t0
                print(f"[{label}] {done}/{n} el={el:.0f}s eta={(n-done)/(done/el):.0f}s", flush=True)
    R_map = results.reshape(len(WLS), len(ANGLES))
    out = os.path.join(outdir, f"map_{label}.npz")
    np.savez(out, wavelength_nm=WLS, angle_deg=ANGLES, R=R_map, nh=NH,
             phi=phi, psi=psi, note=note, computed=True,
             period_nm=P, radius_nm=R, height_nm=H, disc=DISC)
    print(f"[{label}] done {time.time()-t0:.0f}s -> {out}", flush=True)
    return R_map


def write_mirror(label, src_label, R_map, phi, psi, src_note):
    outdir = os.path.join(OUTBASE, label)
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"map_{label}.npz")
    np.savez(out, wavelength_nm=WLS, angle_deg=ANGLES, R=R_map, nh=NH,
             phi=phi, psi=psi, note=f"{src_note}; identical to {src_label} by C4v symmetry",
             computed=False, mirror_of=src_label,
             period_nm=P, radius_nm=R, height_nm=H, disc=DISC)
    print(f"[{label}] written as C4v mirror of {src_label} -> {out}", flush=True)


if __name__ == "__main__":
    os.makedirs(OUTBASE, exist_ok=True)
    print(f"Output base: {OUTBASE}", flush=True)
    print(f"Geometry P={P} R={R} H={H} nm; grid {len(WLS)}x{len(ANGLES)}={len(WLS)*len(ANGLES)} per map", flush=True)
    t0 = time.time()
    maps = {}
    for label, cfg in COMPUTED.items():
        maps[label] = run_one(label, cfg["phi"], cfg["psi"], cfg["note"])
    for label, cfg in MIRRORED.items():
        src = cfg["source"]
        write_mirror(label, src, maps[src], cfg["phi"], cfg["psi"], COMPUTED[src]["note"])
    print(f"ALL DONE in {time.time()-t0:.0f}s", flush=True)
