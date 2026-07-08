# Bare a-Si disk BIC metasurface — angle-resolved reflectance (RCWA / nannos)

Classic RCWA reflectance energy-momentum maps of the **bare** amorphous-silicon
disk metasurface from Marangi et al., *"Enhanced cooperativity of
J-exciton-polaritons in dielectric BIC metasurfaces"*, arXiv:2511.08103
(reproducing the bare-metasurface panel Fig. 2a). **No TDBC / J-aggregate /
molecule coating** — Si disks in air on a fused-silica substrate only.

Generated 2026-07-08. Pipeline: `simulation_scripts/Si_disk_metasurface/`.

## Geometry
- Square lattice, period **P = 320 nm** (both axes).
- a-Si cylinder (disk): **radius R = 80 nm**, **height H = 88 nm**
  (paper ranges: R = 70–90 nm, H = 85–90 nm depending on cavity detuning).
- Stack (top→bottom): air (semi-∞) / patterned a-Si disk layer (88 nm) /
  fused-silica substrate (semi-∞).

## Materials
- **a-Si**: user-measured complex index n + i·k (`material_data/a_Si_complex_refr_index.txt`,
  ellipsometry, ~381–893 nm), linearly interpolated. This film is **lossy** in the
  visible (k ≈ 0.09 at 560 nm, ≈ 0.06 at 590 nm), which caps the material-limited
  Q at ~25–30 — see caveat below.
- **Fused silica**: Malitson 1965 Sellmeier (lossless).
- **Air**: n = 1.

## Simulation
- nannos RCWA, **nh = 140** harmonics, discretization 256×256,
  `formulation="tangent"`. (nh convergence check: the BIC peak R settles by
  nh≈180 — nh 180 vs 220 agree to <0.5 %; nh=140 used here reproduces the BIC
  **position/dispersion exactly** and is ~10–12 % high in absolute peak
  amplitude. Rerun `sweep.py 180 2 1` for a fully amplitude-converged, finer map.)
- Reported quantity: **total reflectance** R (sum over reflected orders). The
  structure is subwavelength across 500–800 nm in both air and SiO₂, so total R
  equals the zeroth-order specular reflectance.
- Grid: λ = 500–800 nm @ 3 nm × θ = −40…+40° @ 1.5° (101 × 54 per map).
- Output-port slit convention: **"ky namely kx=0"** = tilt along y (kx fixed at 0,
  sweep ky); **"kx namely ky=0"** = tilt along x (ky fixed at 0, sweep kx).

## The four maps (subfolders)
| folder | input pol | sweep (slit) | φ, ψ | character |
|--------|-----------|--------------|------|-----------|
| `1a_ypol_ky-sweep_kx0` | y | ky-sweep (kx=0) | 90, 0  | TM-like — BIC **dark** |
| `1b_ypol_kx-sweep_ky0` | y | kx-sweep (ky=0) | 0, 90  | TE-like — BIC **bright** |
| `2a_xpol_ky-sweep_kx0` | x | ky-sweep (kx=0) | 90, 90 | TE-like — BIC **bright** |
| `2b_xpol_kx-sweep_ky0` | x | kx-sweep (ky=0) | 0, 0   | TM-like — BIC **dark** |

Each subfolder holds `map_<label>.npz` (arrays `wavelength_nm`, `angle_deg`, `R`,
plus geometry/settings metadata) and a two-panel `map_<label>.png`
(left: full range; right: contrast-enhanced, vmax=0.06, to reveal the weak BIC).

## Important: C4v symmetry → only 2 distinct maps
An **ideal circular disk** on a square lattice is C4v-symmetric. A 90° rotation
maps x↔y and kx↔ky, so **exactly**:
- `2a` (x, ky-sweep) **≡** `1b` (y, kx-sweep)
- `2b` (x, kx-sweep) **≡** `1a` (y, ky-sweep)

verified numerically to ~1×10⁻¹⁴. The cross-polarization maps therefore carry no
independent information for an ideal disk. Maps `1a` and `1b` were computed
directly; `2a` and `2b` were written from the exact symmetry relation (flagged
`computed=False`, `mirror_of=...` in their npz). A real, slightly asymmetric
(C2) fabricated pillar — as in the paper — would lift this degeneracy and make
the cross-polarization maps genuinely independent.

## Physics notes
- **BIC** (symmetry-protected, out-of-plane MD): appears only in the **TE**
  configuration (E ⟂ plane of incidence → `1b`/`2a`), as a dispersive Fano band
  redshifting with |angle| (~540 nm near normal → ~575 nm at 20°), narrowing and
  fading toward θ=0° where it decouples. It is **weak and broad here** (contrast
  ΔR ~ 0.02, Q ~ tens) because the measured a-Si is absorptive — a faithful
  consequence of the supplied index, not a numerical artifact. The paper's much
  sharper BIC implies a lower-loss a-Si than this ellipsometry table.
- **Mie band**: strong, rises toward ~500–520 nm, and is likewise strongly
  polarization/geometry dependent (bright in TE `1b`/`2a`, weak in TM `1a`/`2b`).
