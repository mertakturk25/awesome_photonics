import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multipole import multipole_at


def _worker(args):
    wl, theta = args
    res = multipole_at(wl, theta, 90.0, 90.0, nh=100, disc=(128, 128), nz=11)
    return res["abs_p2"], res["abs_T2"], res["abs_peff2"]


if __name__ == "__main__":
    jobs = [(wl, th) for wl in [560, 570, 580, 588.5, 600, 610, 620, 630] for th in [0, 10]]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        results = list(ex.map(_worker, jobs))
    dt = time.time() - t0
    print(f"{len(jobs)} points in {dt:.1f}s -> {dt/len(jobs):.2f} s/point effective (8-way)")
