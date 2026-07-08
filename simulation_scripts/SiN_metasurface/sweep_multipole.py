"""
Compute the ED/toroidal multipole decomposition (|p|^2, |T|^2, |p_eff|^2)
across the same ky-sweep (kx=0, Ex-polarized) grid as the reflectance map,
using the reduced-cost settings agreed with the user:
  wavelength: 450-750nm @ 5nm
  angle: -40 to 40 deg @ 2deg
  nz: 5 z-slices through the 148nm SiN layer

Run as: python sweep_multipole.py
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multipole import multipole_at

NH = 100
DISC = (128, 128)
NZ = 5
WL_MIN, WL_MAX, WL_STEP = 450.0, 750.0, 5.0
ANGLE_MIN, ANGLE_MAX, ANGLE_STEP = -40.0, 40.0, 2.0

WLS = np.arange(WL_MIN, WL_MAX + 1e-6, WL_STEP)
ANGLES = np.arange(ANGLE_MIN, ANGLE_MAX + 1e-6, ANGLE_STEP)

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def _worker(args):
    wl, theta = args
    res = multipole_at(wl, theta, 90.0, 90.0, nh=NH, disc=DISC, nz=NZ)
    return res["R"], res["abs_p2"], res["abs_T2"], res["abs_peff2"]


if __name__ == "__main__":
    jobs = [(wl, th) for wl in WLS for th in ANGLES]
    n = len(jobs)
    print(f"{n} points, nh={NH}, nz={NZ}, wl=[{WL_MIN},{WL_MAX}]@{WL_STEP}, angle=[{ANGLE_MIN},{ANGLE_MAX}]@{ANGLE_STEP}", flush=True)
    t0 = time.time()
    R = np.zeros(n)
    P2 = np.zeros(n)
    T2 = np.zeros(n)
    PEFF2 = np.zeros(n)
    done = 0
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        futs = {ex.submit(_worker, j): i for i, j in enumerate(jobs)}
        for fut in as_completed(futs):
            i = futs[fut]
            r, p2, t2, peff2 = fut.result()
            R[i], P2[i], T2[i], PEFF2[i] = r, p2, t2, peff2
            done += 1
            if done % 100 == 0:
                elapsed = time.time() - t0
                rate = done / elapsed
                eta = (n - done) / rate
                print(f"{done}/{n}  elapsed={elapsed:.0f}s  eta={eta:.0f}s", flush=True)
    shape = (len(WLS), len(ANGLES))
    np.savez(
        os.path.join(OUTDIR, "map_ky_multipole.npz"),
        wavelength_nm=WLS,
        angle_deg=ANGLES,
        R=R.reshape(shape),
        abs_p2=P2.reshape(shape),
        abs_T2=T2.reshape(shape),
        abs_peff2=PEFF2.reshape(shape),
        nh=NH,
        nz=NZ,
    )
    print(f"done in {time.time()-t0:.0f}s -> map_ky_multipole.npz", flush=True)
