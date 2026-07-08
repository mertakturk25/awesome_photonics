"""
RCWA (nannos) model of the BARE amorphous-silicon disk BIC metasurface
(Marangi et al., arXiv:2511.08103, Fig. 2a -- no J-aggregate / molecule coating).

Geometry (confirmed with the user):
  P  = 320 nm      (square lattice period, both axes)
  R  = 80 nm       (a-Si disk radius; paper range 70-90 nm vs detuning)
  H  = 88 nm       (a-Si disk height; paper range 85-90 nm)
  Stack (top -> bottom):
    Air (semi-infinite)  /  patterned a-Si layer (disk in air, H=88nm)  /
    fused-silica substrate (semi-infinite)
  Single circular disk centred in the unit cell.

Note on symmetry: a single circular disk on a square lattice is C4v, so the
BIC couples identically to x- and y-polarization at Gamma and the two in-plane
sweep directions are related by a 90-deg rotation. The paper's "C2" label refers
to their fabricated (slightly asymmetric) pillar; here we model an ideal disk as
requested.

Illumination (nannos PlaneWave(angles=(theta,phi,psi)), amplitude
cx = cos(psi)cos(theta)cos(phi) - sin(psi)sin(phi),
cy = cos(psi)cos(theta)sin(phi) + sin(psi)cos(phi)):
  E along y, ky-sweep (kx=0):  phi=90, psi=0   -> in-plane E purely y
  E along y, kx-sweep (ky=0):  phi=0,  psi=90  -> in-plane E purely y
  E along x, ky-sweep (kx=0):  phi=90, psi=90  -> in-plane E purely x
  E along x, kx-sweep (ky=0):  phi=0,  psi=0   -> in-plane E purely x

Reported quantity: total reflectance (sum over reflected diffraction orders).
For P=320 nm the structure is subwavelength across 500-800 nm in both air and
the SiO2 substrate (no propagating diffraction orders), so total R equals the
zeroth-order specular reflectance.
"""

import os
import sys
import numpy as np
import nannos as nn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from materials import asi_n, sio2_n

P = 320.0
R = 80.0
H = 88.0
CENTER = (P / 2, P / 2)

# (phi, psi) for each requested (polarization, sweep) combination
CONFIG = {
    "1a_ypol_ky-sweep_kx0": dict(phi=90.0, psi=0.0),
    "1b_ypol_kx-sweep_ky0": dict(phi=0.0, psi=90.0),
    "2a_xpol_ky-sweep_kx0": dict(phi=90.0, psi=90.0),
    "2b_xpol_kx-sweep_ky0": dict(phi=0.0, psi=0.0),
}


def build_stack(lattice, wl_nm):
    n_asi = complex(asi_n(wl_nm))
    n_sio2 = complex(sio2_n(wl_nm))

    eps_air = 1.0 + 0j
    eps_asi = n_asi**2
    eps_sio2 = n_sio2**2

    superstrate = lattice.Layer("Air", epsilon=eps_air)

    grid = lattice.ones() * eps_air
    disk = lattice.circle(CENTER, R)
    grid[disk] = eps_asi
    patterned = lattice.Layer("aSi", thickness=H)
    patterned.epsilon = grid

    substrate = lattice.Layer("SiO2", epsilon=eps_sio2)
    return [superstrate, patterned, substrate]


def total_R(wl_nm, theta_deg, phi_deg, psi_deg, nh=200, disc=(256, 256), formulation="tangent"):
    lattice = nn.Lattice(((P, 0), (0, P)), disc)
    stack = build_stack(lattice, wl_nm)
    pw = nn.PlaneWave(wavelength=wl_nm, angles=(theta_deg, phi_deg, psi_deg))
    sim = nn.Simulation(stack, pw, nh=nh, formulation=formulation)
    R, _T = sim.diffraction_efficiencies()
    return float(R)


if __name__ == "__main__":
    # quick self-test at normal incidence
    for wl in (560, 590, 620):
        r = total_R(wl, 0.0, 0.0, 0.0, nh=100)
        na = asi_n(wl)
        print(f"wl={wl}nm  n_aSi={na.real:.3f}+{na.imag:.3f}i  R(normal)={r:.4f}")
