"""
Material dispersion for the Crotti et al. (Light: Sci. Appl. 2024) AlGaAs
nanowire metasurface.

Data source: Papatryfonos et al., "Refractive indices of MBE-grown
AlxGa1-xAs ternary alloys in the transparent wavelength region", AIP Adv. 11,
025327 (2021), as hosted in the refractiveindex.info database
(polyanskiy/refractiveindex.info-database, CC0), cloned locally.

The database only ships discrete measured compositions (0%, 9.7%, 21.9%, ...
Al). The paper's wires are Al0.18Ga0.82As (nominal 18%, "estimated 18.5%" in
the text, consistent with the stated bandgap EG = 1.6528 eV via the standard
Eg(x) = 1.424 + 1.247x direct-gap relation -> x = 0.1835). Since 0.185 falls
between the two nearest measured points (9.7% and 21.9%), n and k are
linearly interpolated in composition at each wavelength. This is a standard
approximation for closely-spaced III-V ternary compositions away from sharp
band-edge features; because our window (650-800 nm) straddles the band edge,
this is flagged explicitly to the user as an approximation (see README note
at bottom of this file / final report).
"""

import os
import numpy as np
import yaml

RII_DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "material_data")


def _load_papatryfonos(fname):
    with open(f"{RII_DB}\\{fname}", encoding="utf-8") as f:
        d = yaml.safe_load(f)
    lines = d["DATA"][0]["data"].strip().split("\n")
    arr = np.array([[float(v) for v in l.split()] for l in lines])
    wl_um, n, k = arr[:, 0], arr[:, 1], arr[:, 2]
    return wl_um * 1000.0, n, k  # nm, n, k


def algaas_nk(wl_nm, x=0.185):
    """Interpolated Al_x Ga_(1-x) As complex refractive index N = n + ik."""
    wl9, n9, k9 = _load_papatryfonos("Papatryfonos-9.7.yml")
    wl22, n22, k22 = _load_papatryfonos("Papatryfonos-21.9.yml")
    x9, x22 = 0.097, 0.219
    t = (x - x9) / (x22 - x9)

    n9i = np.interp(wl_nm, wl9, n9)
    k9i = np.interp(wl_nm, wl9, k9)
    n22i = np.interp(wl_nm, wl22, n22)
    k22i = np.interp(wl_nm, wl22, k22)

    n = (1 - t) * n9i + t * n22i
    k = (1 - t) * k9i + t * k22i
    return n + 1j * k


def gaas_nk(wl_nm):
    """Pure GaAs (x=0) substrate, same Papatryfonos measurement series."""
    wl0, n0, k0 = _load_papatryfonos("Papatryfonos-0.yml")
    n = np.interp(wl_nm, wl0, n0)
    k = np.interp(wl_nm, wl0, k0)
    return n + 1j * k


def alox_n(wl_nm):
    """AlOx buffer: constant, lossless n=1.6 per paper's SI (S1)."""
    return np.full(np.shape(wl_nm), 1.6 + 0j)


if __name__ == "__main__":
    wl = np.arange(650, 801, 5.0)
    N_algaas = algaas_nk(wl)
    N_gaas = gaas_nk(wl)
    Eg_formula = 1.424 + 1.247 * 0.185
    lam_g_formula = 1239.84 / Eg_formula
    print(f"Target composition x=0.185, formula Eg={Eg_formula:.4f} eV -> lambda_g={lam_g_formula:.1f} nm")
    print(f"Paper states EG = 1.6528 eV -> lambda_g = {1239.84/1.6528:.1f} nm")
    print()
    print(" wl(nm)   n_AlGaAs   k_AlGaAs   n_GaAs   k_GaAs")
    for w, Na, Ng in zip(wl, N_algaas, N_gaas):
        print(f" {w:6.1f}  {Na.real:8.4f}  {Na.imag:8.4f}  {Ng.real:7.4f}  {Ng.imag:7.4f}")
