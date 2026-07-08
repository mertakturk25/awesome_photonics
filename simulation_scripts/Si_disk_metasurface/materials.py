"""
Material dispersion for the bare amorphous-silicon (a-Si) disk metasurface
(Marangi et al., "Enhanced cooperativity of J-exciton-polaritons in dielectric
BIC metasurfaces", arXiv:2511.08103 -- bare BIC metasurface of Fig. 2a, WITHOUT
the TDBC J-aggregate coating).

  - a-Si disks: complex refractive index n + i*k from the user's own measured
    ellipsometry table `material_data/a_Si_complex_refr_index.txt`
    (columns: wavelength[nm]  n  k), covering ~381-893 nm. Linearly
    interpolated in n and k separately. This is lossy (k>0) across the
    visible, which sets the finite Q of the BIC.
  - Fused-silica (quartz) substrate: Malitson 1965 Sellmeier, lossless.
  - Air superstrate: n = 1.

The a-Si table was measured on a 90-nm a-Si film (Si_90nm_300C ellipsometry) and
confirmed to match the "new" a_Si_complex_refr_index.txt exactly.
"""

import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ASI_FILE = os.path.join(_HERE, "material_data", "a_Si_complex_refr_index.txt")


def _load_asi():
    wl, n, k = [], [], []
    with open(_ASI_FILE) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                w, nn, kk = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue  # skip header line
            wl.append(w)
            n.append(nn)
            k.append(kk)
    idx = np.argsort(wl)
    return np.array(wl)[idx], np.array(n)[idx], np.array(k)[idx]


_ASI_WL, _ASI_N, _ASI_K = _load_asi()


def asi_n(wl_nm):
    """Complex a-Si refractive index n + i*k at wl_nm (nm), tabulated/interpolated."""
    wl = np.asarray(wl_nm, dtype=float)
    n = np.interp(wl, _ASI_WL, _ASI_N)
    k = np.interp(wl, _ASI_WL, _ASI_K)
    return n + 1j * k


def _sellmeier(wl_nm, coeffs):
    lam_um = np.asarray(wl_nm, dtype=float) / 1000.0
    n2 = 1.0 + coeffs[0]
    for i in range(1, len(coeffs), 2):
        ci, gi = coeffs[i], coeffs[i + 1]
        n2 = n2 + ci * lam_um**2 / (lam_um**2 - gi**2)
    return np.sqrt(n2)


# Malitson 1965, J. Opt. Soc. Am. 55, 1205 -- fused silica, valid 0.21-6.7 um
_SIO2_COEFFS = [0, 0.6961663, 0.0684043, 0.4079426, 0.1162414, 0.8974794, 9.896161]


def sio2_n(wl_nm):
    """Fused-silica refractive index (real, lossless) at wl_nm (nm)."""
    return _sellmeier(wl_nm, _SIO2_COEFFS) + 0j


if __name__ == "__main__":
    print(f"a-Si table: {_ASI_WL.min():.1f}-{_ASI_WL.max():.1f} nm, {len(_ASI_WL)} points")
    for w in np.arange(500, 801, 50.0):
        nasi = asi_n(w)
        print(f"{w:6.1f} nm   n_aSi={nasi.real:.4f}+{nasi.imag:.4f}i   n_SiO2={sio2_n(w).real:.4f}")
