# a-Si disk metasurface — disk shrunk to R=70 nm to separate BIC from Mie (high-res)

Variant of the bare a-Si disk metasurface (Marangi et al., arXiv:2511.08103, no
coating) that separates the BIC from the Mie resonance by **shrinking the disk**
(keeping the period fixed), computed at **higher resolution**. This is the cleaner
of the two separation strategies explored (the other, enlarging the period to
P=380, is in `..._P380-BIC-Mie-split/`).

## Geometry
- Square lattice **P = 320 nm** (unchanged from baseline).
- a-Si disk **R = 70 nm** (baseline 80 nm), **H = 88 nm**.
- air / a-Si disk in air / fused-silica substrate.

## Why shrink the disk (vs enlarging the period)
- **Mie** is a single-disk mode set by R, H → shrinking the disk **blueshifts** it
  (R=80→Mie~488 nm, R=70→~475 nm) and weakens it.
- **BIC** is a lattice mode (~n_eff·P) → fixing P keeps it at ~530–540 nm.
- **Substrate Rayleigh/Wood anomaly** at n_SiO2·P ≈ **468 nm** → fixing P keeps it
  DOWN near the Mie, **out of the BIC region**.

So the Mie moves down toward ~475 nm (next to the 468 nm Rayleigh line) while the
BIC stays ~530 nm and disperses upward — the BIC band is cleanly separated above,
with no diffraction clutter in its spectral region (the key advantage over the
P=380 route, where the Rayleigh anomaly lands right at the BIC's Γ point).
R=70 was chosen from a scan (R=80/70/62/55) as the best separation-vs-visibility
tradeoff: smaller R separates more but the Mie fades fast (R=62→0.13, R=55→0.03).

## Feature layout
- **Mie**: bright band ~470–500 nm (brightest near normal), just above the Rayleigh
  line at ~468 nm.
- **BIC** (TE only, maps 1b/2a): dispersive band from ~530 nm at Γ, redshifting with
  |angle|; absorption-broadened/weak (view the grayscale/enhanced PNGs).
- Maps 1a/2b (TM) show the Mie/Rayleigh but the BIC is dark.

## Resolution / settings
- Grid **234 × 81** per map: λ = 400–750 nm @ **1.5 nm** × θ = −40…+40° @ **1°**
  (~4× the pixels of the baseline/P380 runs).
- nannos RCWA, nh=140, disc 256², formulation="tangent".
- Same 4-map layout and C4v-symmetry handling as the baseline (2a≡1b, 2b≡1a via the
  exact symmetry relation; see the baseline folder's README for the φ/ψ table).

Each subfolder holds map_*.npz, a two-panel enhanced PNG, and a single-panel
`_gray.png`; plus overview_4maps_enhanced.png and this README.
