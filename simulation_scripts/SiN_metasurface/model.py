"""
RCWA (nannos) model of the Si3N4 nanohole metasurface (PhD thesis Ch.4.3):
  px = py = 380 nm  (square lattice)
  D  = 150 nm       (hole diameter)
  hz = 148 nm       (film thickness, holes etched fully through)
  delta_x = 200 nm, delta_y = 100 nm  (centre-to-centre offset between the
                                       two holes in the unit cell)
  Stack (top->bottom): Air (semi-inf) / patterned SiN layer (148nm) /
                        SiO2 fused silica (semi-inf)

Hole placement: the two holes are placed symmetrically about the unit-cell
centre, offset by +-(delta_x/2, delta_y/2), i.e.
  hole A center = (px/2 - delta_x/2, py/2 - delta_y/2)
  hole B center = (px/2 + delta_x/2, py/2 + delta_y/2)
This keeps both holes fully inside the cell with no periodic wraparound and
reproduces the diagonal centre-to-centre vector (delta_x, delta_y) shown in
the thesis figure.

Illumination: E field fixed along x (thesis: resonance excited for
E parallel to x, vanishes along the dimer diagonal). For a 2D lattice, using
nannos PlaneWave(angles=(theta,phi,psi)):
  phi=90 -> wavevector tilts in the y-z plane (kx=0, "ky-sweep")
  psi=90 -> cx=-sin(phi), cy=cos(phi)*cos(theta)... => at phi=90 this gives
            pure E along x for all theta (same algebra as the AlGaAs project's
            ky-sweep TM convention, which happens to coincide with pure
            Ex-polarization here too).
"""

import os
import sys
import numpy as np
import nannos as nn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from materials import sin_n, sio2_n

PX = 380.0
PY = 380.0
D = 150.0
HZ = 148.0
DX = 200.0
DY = 100.0

HOLE_A = (PX / 2 - DX / 2, PY / 2 - DY / 2)
HOLE_B = (PX / 2 + DX / 2, PY / 2 + DY / 2)
RADIUS = D / 2


def build_stack(lattice, wl_nm):
    n_sin = complex(sin_n(np.array([wl_nm]))[0])
    n_sio2 = complex(sio2_n(np.array([wl_nm]))[0])

    eps_air = 1.0 + 0j
    eps_sin = n_sin**2
    eps_sio2 = n_sio2**2

    superstrate = lattice.Layer("Air", epsilon=eps_air)

    grid = lattice.ones() * eps_sin
    holeA_mask = lattice.circle(HOLE_A, RADIUS)
    holeB_mask = lattice.circle(HOLE_B, RADIUS)
    grid[holeA_mask] = eps_air
    grid[holeB_mask] = eps_air
    patterned = lattice.Layer("SiN", thickness=HZ)
    patterned.epsilon = grid

    substrate = lattice.Layer("SiO2", epsilon=eps_sio2)

    return [superstrate, patterned, substrate]


def specular_R(wl_nm, theta_deg, phi_deg, psi_deg, nh=200, disc=(128, 128), formulation="tangent"):
    lattice = nn.Lattice(((PX, 0), (0, PY)), disc)
    stack = build_stack(lattice, wl_nm)
    pw = nn.PlaneWave(wavelength=wl_nm, angles=(theta_deg, phi_deg, psi_deg))
    sim = nn.Simulation(stack, pw, nh=nh, formulation=formulation)
    Ri, _Ti = sim.diffraction_efficiencies(orders=True)
    R0 = sim.get_order(Ri, (0, 0))
    return float(R0)


if __name__ == "__main__":
    print(f"hole A center = {HOLE_A}, hole B center = {HOLE_B}, radius={RADIUS}")
    dist = np.hypot(HOLE_B[0] - HOLE_A[0], HOLE_B[1] - HOLE_A[1])
    print(f"centre-to-centre distance = {dist:.1f} nm (2*radius = {2*RADIUS} nm, gap = {dist-2*RADIUS:.1f} nm)")
