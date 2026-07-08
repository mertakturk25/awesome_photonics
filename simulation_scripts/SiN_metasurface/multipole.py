"""
Electric dipole (ED) / toroidal dipole (TD) multipole decomposition of the
induced polarization current inside the air holes (the "meta-atom" = absence
of SiN), to look for anapole-like ED-TD cancellation.

Method: run the same RCWA simulation as model.py, extract the E-field on the
2D lattice grid at several z-slices through the 148nm SiN layer, mask to the
hole region (air, where the induced polarization current relative to the SiN
host is nonzero), and integrate the standard multipole moment formulas over
that volume (origin = unit-cell center = centroid of the two-hole dimer).

Since the induced current everywhere in the hole is J(r) = C * E(r) with a
SINGLE complex scalar C = -i*omega*eps0*(eps_air-eps_SiN) (Delta-eps is
spatially uniform inside the hole), C factors out of the anapole condition
p + i*k0*T = 0 entirely -- so we work with the "reduced" moments
  p' = (1/(i*omega)) * Integral[E] dV
  T' = (1/(10c)) * Integral[(r.E) r - 2 r^2 E] dV
and check p' + i*k0*T' for cancellation (C cancels out of this sum's zero).
"""

import os
import sys
import numpy as np
import nannos as nn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from materials import sin_n, sio2_n
from model import PX, PY, D, HZ, HOLE_A, HOLE_B, RADIUS

C_LIGHT = 299792458e9  # nm/s


def multipole_at(wl_nm, theta_deg, phi_deg, psi_deg, nh=100, disc=(128, 128), nz=11):
    n_sin = complex(sin_n(np.array([wl_nm]))[0])
    n_sio2 = complex(sio2_n(np.array([wl_nm]))[0])
    eps_air = 1.0 + 0j
    eps_sin = n_sin**2
    eps_sio2 = n_sio2**2

    lattice = nn.Lattice(((PX, 0), (0, PY)), disc)
    superstrate = lattice.Layer("Air", epsilon=eps_air)
    grid = lattice.ones() * eps_sin
    holeA_mask = lattice.circle(HOLE_A, RADIUS)
    holeB_mask = lattice.circle(HOLE_B, RADIUS)
    grid[holeA_mask] = eps_air
    grid[holeB_mask] = eps_air
    patterned = lattice.Layer("SiN", thickness=HZ)
    patterned.epsilon = grid
    substrate = lattice.Layer("SiO2", epsilon=eps_sio2)
    stack = [superstrate, patterned, substrate]

    pw = nn.PlaneWave(wavelength=wl_nm, angles=(theta_deg, phi_deg, psi_deg))
    sim = nn.Simulation(stack, pw, nh=nh, formulation="tangent")
    R0 = sim.get_order(sim.diffraction_efficiencies(orders=True)[0], (0, 0))

    hole_mask = holeA_mask | holeB_mask
    ny, nx = hole_mask.shape
    xs = np.linspace(0, PX, nx, endpoint=False)
    ys = np.linspace(0, PY, ny, endpoint=False)
    X, Y = np.meshgrid(xs, ys)
    ox, oy = PX / 2, PY / 2  # unit-cell center = dimer centroid
    Xr, Yr = X - ox, Y - oy

    zs = np.linspace(0, HZ, nz)
    dz = zs[1] - zs[0]
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    omega = 2 * np.pi * C_LIGHT / wl_nm
    k0 = 2 * np.pi / wl_nm

    p_int = np.zeros(3, dtype=complex)
    T_int = np.zeros(3, dtype=complex)

    for z in zs:
        E = sim.get_Efield_grid(1, z=z, component="all")  # layer index 1 = SiN
        Ex, Ey, Ez = E[0], E[1], E[2]
        m = hole_mask
        Exm, Eym, Ezm = Ex[m], Ey[m], Ez[m]
        Xm, Ym = Xr[m], Yr[m]
        Zm = np.full_like(Xm, z - HZ / 2)  # z measured from mid-plane of layer

        # electric dipole integrand: just E
        p_int[0] += Exm.sum() * dx * dy * dz
        p_int[1] += Eym.sum() * dx * dy * dz
        p_int[2] += Ezm.sum() * dx * dy * dz

        # toroidal dipole integrand: (r.E) r - 2 r^2 E
        rE = Xm * Exm + Ym * Eym + Zm * Ezm
        r2 = Xm**2 + Ym**2 + Zm**2
        Tx = rE * Xm - 2 * r2 * Exm
        Ty = rE * Ym - 2 * r2 * Eym
        Tz = rE * Zm - 2 * r2 * Ezm
        T_int[0] += Tx.sum() * dx * dy * dz
        T_int[1] += Ty.sum() * dx * dy * dz
        T_int[2] += Tz.sum() * dx * dy * dz

    p_prime = p_int / (1j * omega)
    T_prime = T_int / (10 * C_LIGHT)

    p_eff = p_prime + 1j * k0 * T_prime

    return {
        "R": float(R0),
        "p": p_prime,
        "T": T_prime,
        "p_eff": p_eff,
        "abs_p2": float(np.sum(np.abs(p_prime) ** 2)),
        "abs_T2": float(np.sum(np.abs(T_prime) ** 2)),
        "abs_peff2": float(np.sum(np.abs(p_eff) ** 2)),
    }


if __name__ == "__main__":
    import time

    for wl in [560, 575, 585, 588.5, 590, 595, 600]:
        t0 = time.time()
        res = multipole_at(wl, 0.0, 90.0, 90.0, nh=100, disc=(128, 128), nz=11)
        dt = time.time() - t0
        ratio = res["abs_peff2"] / max(res["abs_p2"], 1e-30)
        print(
            f"wl={wl:6.1f}  R={res['R']:.4f}  |p|^2={res['abs_p2']:.3e}  "
            f"|T|^2={res['abs_T2']:.3e}  |p_eff|^2={res['abs_peff2']:.3e}  "
            f"|p_eff|^2/|p|^2={ratio:.4f}  time={dt:.2f}s"
        )
