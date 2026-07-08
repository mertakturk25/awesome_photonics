"""
RCWA (nannos) model of the Crotti et al. AlGaAs nanowire metasurface,
fitted-geometry configuration (Methods + SI S4.1):
  P = 400 nm (period, grating axis = x, perpendicular to wires)
  W = 150 nm (wire width)
  H = 400 nm (wire height)
  70 nm fillet radius rounding the top corners of the wire
  AlOx buffer = 890 nm, n = 1.6 (lossless, per SI S1)
  GaAs semi-infinite substrate
  Air semi-infinite superstrate

Wires are invariant along y (nannos' non-periodic in-plane direction for a
mono-periodic lattice). TM polarization = E perpendicular to the wires
(no E_y component), matching the paper's convention (Fig. 1a).

For a mono-periodic nannos Lattice, PlaneWave(angles=(theta, phi, psi)) with
theta = polar angle, phi = azimuthal angle, psi = polarization angle gives a
field amplitude direction (cx,cy,cz) with
  cy = cos(psi)*cos(theta)*sin(phi) + sin(psi)*cos(phi)
Setting cy = 0 for all theta requires:
  phi = 0   -> psi = 0   (this is the plane perpendicular to the wires;
                           matches SI Fig. S2's existing 1D dispersion sweep)
  phi = 90  -> psi = 90  (conical mount, tilting along the wire axis)
This was verified against nannos' own 1D grating example (psi=0 <-> "TM",
psi=90 <-> "TE" at phi=0).
"""

import os
import sys
import numpy as np
import nannos as nn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from materials import algaas_nk, gaas_nk, alox_n

P = 400.0
W = 150.0
H = 400.0
FILLET_R = 70.0
ALOX_T = 890.0
N_FILLET_SLABS = 6


def fillet_width(s):
    """Width of the wire at distance s below the very top (s in [0, FILLET_R])."""
    half = (W / 2 - FILLET_R) + np.sqrt(FILLET_R**2 - (FILLET_R - s) ** 2)
    return 2 * half


def wire_slabs():
    """Return list of (thickness, width) from TOP of wire to BOTTOM."""
    slabs = []
    ds = FILLET_R / N_FILLET_SLABS
    for i in range(N_FILLET_SLABS):
        s_mid = (i + 0.5) * ds  # distance below top, at slab midpoint
        slabs.append((ds, fillet_width(s_mid)))
    bulk_thickness = H - FILLET_R
    slabs.append((bulk_thickness, W))
    return slabs


def build_stack(lattice, wl_nm, nh_disc=512):
    n_algaas = complex(algaas_nk(np.array([wl_nm]))[0])
    n_gaas = complex(gaas_nk(np.array([wl_nm]))[0])
    n_alox = 1.6 + 0j

    eps_air = 1.0 + 0j
    eps_algaas = n_algaas**2
    eps_gaas = n_gaas**2
    eps_alox = n_alox**2

    superstrate = lattice.Layer("Air", epsilon=eps_air)

    wire_layers = []
    for idx, (thickness, width) in enumerate(wire_slabs()):
        grid = lattice.ones() * eps_air
        mask = lattice.stripe(P / 2, width)
        grid[mask] = eps_algaas
        layer = lattice.Layer(f"wire_{idx}", thickness=thickness)
        layer.epsilon = grid
        wire_layers.append(layer)

    buffer_layer = lattice.Layer("AlOx", thickness=ALOX_T, epsilon=eps_alox)
    substrate = lattice.Layer("GaAs", epsilon=eps_gaas)

    return [superstrate] + wire_layers + [buffer_layer, substrate]


def specular_R(wl_nm, theta_deg, phi_deg, psi_deg, nh=100, disc=256, formulation="tangent"):
    lattice = nn.Lattice(P, disc)
    stack = build_stack(lattice, wl_nm)
    pw = nn.PlaneWave(wavelength=wl_nm, angles=(theta_deg, phi_deg, psi_deg))
    sim = nn.Simulation(stack, pw, nh=nh, formulation=formulation)
    Ri, _Ti = sim.diffraction_efficiencies(orders=True)
    R0 = sim.get_order(Ri, 0)
    return float(R0)


if __name__ == "__main__":
    # quick sanity print of the fillet staircase
    for t, w in wire_slabs():
        print(f"thickness={t:6.2f} nm   width={w:6.2f} nm")
