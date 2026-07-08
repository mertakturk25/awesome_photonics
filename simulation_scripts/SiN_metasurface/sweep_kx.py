"""
Compute the Ex-polarized static reflectance energy-momentum map for the
Si3N4 nanohole metasurface, ky-sweep (kx=0): phi=90, psi=90, theta in
[-60,60] deg, over the wavelength range given on the command line.

Run as: python sweep.py [nh] [wl_step] [angle_step] [wl_min] [wl_max]
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

NH = int(sys.argv[1]) if len(sys.argv) > 1 else 100
WL_STEP = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
ANGLE_STEP = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
WL_MIN = float(sys.argv[4]) if len(sys.argv) > 4 else 400.0
WL_MAX = float(sys.argv[5]) if len(sys.argv) > 5 else 750.0
DISC = (128, 128)

WLS = np.arange(WL_MIN, WL_MAX + 1e-6, WL_STEP)
ANGLES = np.arange(-60.0, 60.0 + 1e-6, ANGLE_STEP)

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def _worker(args):
    wl, theta, phi, psi = args
    return specular_R(wl, theta, phi, psi, nh=NH, disc=DISC)


def run_sweep(phi, psi, label):
    jobs = [(wl, th, phi, psi) for wl in WLS for th in ANGLES]
    n = len(jobs)
    print(f"[{label}] {n} points, nh={NH}, wl=[{WL_MIN},{WL_MAX}]@{WL_STEP}, angle_step={ANGLE_STEP}", flush=True)
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
                print(f"[{label}] {done}/{n}  elapsed={elapsed:.0f}s  eta={eta:.0f}s", flush=True)
    R_map = results.reshape(len(WLS), len(ANGLES))
    np.savez(os.path.join(OUTDIR, f"map_{label}.npz"), wavelength_nm=WLS, angle_deg=ANGLES, R=R_map, nh=NH)
    print(f"[{label}] done in {time.time()-t0:.0f}s -> map_{label}.npz", flush=True)
    return R_map


if __name__ == "__main__":
    print("Starting kx sweep (phi=0, psi=0, ky=0, Ex-polarized)...", flush=True)
    run_sweep(phi=0.0, psi=0.0, label="kx")
    print("ALL DONE", flush=True)
