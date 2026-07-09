# a-Si disk BIC metasurface — period increased to P=380 nm to separate BIC from Mie

Variant of the bare a-Si disk metasurface (Marangi et al., arXiv:2511.08103,
no coating) with the **period increased to P=380 nm** to push the BIC away from
the Mie resonance. Same disk size and materials as the baseline
`Si_disk_20260708_bare-reflectance/`, wavelength window **400–750 nm**.

## Geometry
- Square lattice **P = 380 nm** (baseline was 320 nm).
- a-Si disk **R = 80 nm, H = 88 nm** (unchanged).
- air / a-Si disk in air / fused-silica substrate.

## Why P was increased
The BIC is a **lattice mode** (its wavelength scales ~ n_eff·P), while the Mie is
a **single-disk particle mode** (set by R, H, nearly period-independent). So
increasing P redshifts the BIC while leaving the Mie roughly in place, opening the
spectral gap between them. A parametric scan (P = 320/360/380/400/440) confirmed
this; P=380 is a sweet spot — clear separation while the resonances stay visible.
P≥440 was rejected: the resonances collapse (R<0.1) and the map fills with
diffraction artifacts.

## Feature identification (from TE vs TM spectra, θ=1° and 25°)
- **BIC** (TE only, i.e. maps 1b/2a): dispersive band from ~555 nm at Γ, redshifting
  to ~644 nm at 25° and ~700 nm by 40°. Absent in TM — confirms it is the
  symmetry-protected mode. Weak (R~0.02–0.05, absorption-broadened by the lossy
  a-Si), so view the grayscale/enhanced PNGs.
- **Mie** (seen in TM, maps 1a/2b): ~480–530 nm, only weakly dispersing.
- **Substrate Rayleigh/Wood anomaly** at n_SiO2·P ≈ **555 nm** (flat line in BOTH
  polarizations): a diffraction threshold, NOT a resonance. Because it also scales
  with P, it now sits in-window at the BIC's Γ point. It is the main cosmetic cost
  of enlarging the period.

## Result vs baseline
- Baseline P=320 (window 500–800 nm): BIC_Γ ~540 nm, Mie ~500 nm — nearly merged.
- Here P=380: BIC_Γ ~555 nm and dispersing up to ~700 nm, Mie ~490 nm — the BIC
  band is clearly separated above the Mie across most of the map.

## Files
Same 4-map layout and conventions as the baseline (see that folder's README for
the nannos φ/ψ table and the C4v-symmetry argument: 2a≡1b, 2b≡1a computed via the
exact symmetry relation). Settings: nh=140, disc 256², formulation="tangent",
grid λ 400–750 nm @3 nm × θ ±40° @1.5°. Each subfolder holds map_*.npz, a
two-panel enhanced PNG, and a single-panel `_gray.png`; plus overview + this README.
