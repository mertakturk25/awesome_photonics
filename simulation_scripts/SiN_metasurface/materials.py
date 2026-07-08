"""
Material dispersion for the Si3N4 nanohole metasurface (PhD thesis Ch.4.3):
  - Si3N4 (SiN) film: Luke et al. 2015 Sellmeier ("formula 1" on
    refractiveindex.info), stoichiometric low-loss LPCVD-type film,
    n only (lossless, matching thesis's own description of negligible
    losses / no Fabry-Perot damping across the visible range).
  - Fused silica substrate: Malitson 1965 Sellmeier, the standard
    reference for fused silica, also lossless in this range.
  - Air: n=1.

Both are "formula 1" (Sellmeier) entries: coefficients [c0, c1, c2, c3, c4, ...]
with wavelength in micrometers:
  n(lambda)^2 = 1 + c0 + sum_i c_(2i-1) * lambda^2 / (lambda^2 - c_(2i)^2)
"""

import numpy as np


def _sellmeier(wl_nm, coeffs):
    lam_um = np.asarray(wl_nm, dtype=float) / 1000.0
    c0 = coeffs[0]
    n2 = 1.0 + c0
    for i in range(1, len(coeffs), 2):
        ci, gi = coeffs[i], coeffs[i + 1]
        n2 = n2 + ci * lam_um**2 / (lam_um**2 - gi**2)
    return np.sqrt(n2)


# Luke et al. 2015, Opt. Lett. 40, 4823 -- Si3N4, valid 0.310-5.504 um
_SIN_COEFFS = [0, 3.0249, 0.1353406, 40314, 1239.842]

# Malitson 1965, J. Opt. Soc. Am. 55, 1205 -- fused silica, valid 0.21-6.7 um
_SIO2_COEFFS = [0, 0.6961663, 0.0684043, 0.4079426, 0.1162414, 0.8974794, 9.896161]


def sin_n(wl_nm):
    """Si3N4 refractive index (real, lossless) at wl_nm (nm)."""
    return _sellmeier(wl_nm, _SIN_COEFFS) + 0j


def sio2_n(wl_nm):
    """Fused silica refractive index (real, lossless) at wl_nm (nm)."""
    return _sellmeier(wl_nm, _SIO2_COEFFS) + 0j


if __name__ == "__main__":
    wl = np.arange(400, 751, 25.0)
    for w in wl:
        print(f"{w:6.1f} nm   n_SiN={sin_n(w).real:.4f}   n_SiO2={sio2_n(w).real:.4f}")
