"""
Plain RCWA far-field (specular reflectance) sweep for the SiN nanohole
metasurface, Ex-polarized, either kx-sweep (phi=0,psi=0,ky=0) or ky-sweep
(phi=90,psi=90,kx=0).

Run as: python sweep_farfield.py <kx|ky> <outdir> [nh] [wl_step] [angle_step] [wl_min] [wl_max] [angle_min] [angle_max]
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import specular_R

DIRECTION = sys.argv[1]
OUTDIR = sys.argv[2]
NH = int(sys.argv[3]) if len(sys.argv) > 3 else 100
WL_STEP = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
ANGLE_STEP = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
WL_MIN = float(sys.argv[6]) if len(sys.argv) > 6 else 450.0
WL_MAX = float(sys.argv[7]) if len(sys.argv) > 7 else 750.0
ANGLE_MIN = float(sys.argv[8]) if len(sys.argv) > 8 else -40.0
ANGLE_MAX = float(sys.argv[9]) if len(sys.argv) > 9 else 40.0
DISC = (128, 128)

if DIRECTION == "kx":
    PHI, PSI = 0.0, 0.0
elif DIRECTION == "ky":
    PHI, PSI = 90.0, 90.0
else:
    raise ValueError("direction must be kx or ky")

WLS = np.arange(WL_MIN, WL_MAX + 1e-6, WL_STEP)
ANGLES = np.arange(ANGLE_MIN, ANGLE_MAX + 1e-6, ANGLE_STEP)


def _worker(args):
    wl, theta = args
    return specular_R(wl, theta, PHI, PSI, nh=NH, disc=DISC)


if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)
    jobs = [(wl, th) for wl in WLS for th in ANGLES]
    n = len(jobs)
    print(
        f"[{DIRECTION}] {n} points, nh={NH}, wl=[{WL_MIN},{WL_MAX}]@{WL_STEP}, "
        f"angle=[{ANGLE_MIN},{ANGLE_MAX}]@{ANGLE_STEP}",
        flush=True,
    )
    t0 = time.time()
    results = np.zeros(n)
    done = 0
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        futs = {ex.submit(_worker, j): i for i, j in enumerate(jobs)}
        for fut in as_completed(futs):
            i = futs[fut]
            results[i] = fut.result()
            done += 1
            if done % 500 == 0:
                elapsed = time.time() - t0
                rate = done / elapsed
                eta = (n - done) / rate
                print(f"[{DIRECTION}] {done}/{n}  elapsed={elapsed:.0f}s  eta={eta:.0f}s", flush=True)
    R_map = results.reshape(len(WLS), len(ANGLES))
    outfile = os.path.join(OUTDIR, f"map_{DIRECTION}_farfield.npz")
    np.savez(outfile, wavelength_nm=WLS, angle_deg=ANGLES, R=R_map, nh=NH)
    print(f"[{DIRECTION}] done in {time.time()-t0:.0f}s -> {outfile}", flush=True)
