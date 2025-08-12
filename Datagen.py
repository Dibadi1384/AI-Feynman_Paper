# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Generate data sets (original & π versions) for five physics relations.
# Author: <your-name>
# """

import numpy as np
import os
from math import pi, sqrt

# # ---------------------------------------------------------------------
# # constants
# EPS0 = 8.854_187_817e-12       # F m-1
# MU0  = 4 * pi * 1e-7           # N A-2
# C    = 1 / sqrt(EPS0 * MU0)    # m s-1
# E_CH = 1.602_176_634e-19       # C
# M_E  = 9.109_383_56e-31        # kg
# HBAR = 1.054_571_817e-34       # J s
# EV2J = E_CH                    # 1 eV in joules
# VP   = 343.0                   # m s-1 (speed of sound in air)

# physical constants (fixed)
ME   = 9.10938356e-31      # electron mass [kg]
HBAR = 1.054571817e-34     # reduced Planck constant [J·s]
EV2J = 1.602176634e-19     # 1 eV in joules

EPS0 = 8.854187817e-12     # vacuum permittivity [F/m]
MU0  = 4 * np.pi * 1e-7    # vacuum permeability [N/A²]

RS = np.random.default_rng(42)  # reproducible seed


def save_data(name: str, data: np.ndarray, expr: str):
    """Utility: save data and the symbolic expression."""
    np.savetxt(f"{name}.txt", data, fmt="%.6e", delimiter=" ")
    with open(f"{name}_equation.txt", "w", encoding="utf-8") as fh:
        fh.write(expr + "\n")


# # Field‐energy density  --------------------------------------------
# def field_energy_original(n=10000):
#     E = RS.uniform(0, 5, n)          # V m-1
#     B = RS.uniform(0, 5, n)          # T
#     u = 0.5 * EPS0 * E**2 + 0.5 * B**2 / MU0
#     # constant columns:
#     eps0_col = EPS0 * np.ones(n)
#     mu0_col  = MU0  * np.ones(n)
#     # ensure u is last column
#     data = np.column_stack([E, B, eps0_col, mu0_col, u])
#     expr = "u = 0.5*ε0*E**2 + 0.5*B**2/μ0"
#     save_data("field_energy_density", data, expr)

# def field_energy_pi(n=10):   # <-------------this one worked!
#     pi2 = RS.uniform(0, 10, n)
#     pi1 = 0.5 * (1.0 + pi2**2)
#     data = np.column_stack([pi2, pi1])
#     expr = "π1 = (1 + π2**2) / 2"
#     save_data("field_energy_density_pi", data, expr)


# # Fermi‐gas density of states  -------------------------------------
# def fermi_dos_original(n=10000):
#     V   = RS.uniform(0, 5, n)    # m3
#     Eev = RS.uniform(0, 5, n)    # eV
#     E   = Eev * EV2J                  # J
#     g   = (V / (2 * pi**2)) * (2 * M_E / HBAR**2)**1.5 * np.sqrt(E)
#     # constant columns:   #for the our purpose im gonna asign them a variable
#     me_col   = M_E  * np.ones(n)
#     hbar_col = HBAR * np.ones(n)
#     ev2j_col = EV2J * np.ones(n)
#     # dependent g last
#     data = np.column_stack([V, E, me_col, hbar_col, ev2j_col, g])
#     expr = "g = V/(2π²) * (2me/ħ²)^{3/2} * E^{1/2}"
#     save_data("fermi_gas_DOS", data, expr)

# def fermi_dos_pi(n=10):
#     pi2 = RS.uniform(0.1, 10.0, n)
#     pi1 = sqrt(2) / pi**2 * pi2
#     data = np.column_stack([pi2, pi1])
#     expr = "π1 = √2 / π² * π2"
#     save_data("fermi_gas_DOS_pi", data, expr)


# # Doppler shift  -------------------------------------
# def doppler_original(n=100):
#     u      = RS.uniform(0, 10.0, n)   # m s-1
#     theta  = RS.uniform(0.0, pi, n)   # rad
#     vratio = 1.0 - u / VP * np.cos(theta)
#     # constant column
#     vp_col = VP * np.ones(n)
#     # dependent vratio last
#     data   = np.column_stack([u, theta, vp_col, vratio])
#     expr   = "v_ratio = 1 - u/VP * cosθ"
#     save_data("doppler_shift", data, expr)


# def doppler_pi(n=10):
#     pi2   = RS.uniform(0, 10, n)
#     theta = RS.uniform(0.0, pi, n)
#     pi1   = 1.0 - pi2 * np.cos(theta)
#     data  = np.column_stack([pi2, theta, pi1])
#     expr  = "π1 = 1 - π2 * cosθ"
#     save_data("doppler_shift_pi", data, expr)

# def stokes_pi_groups_data(n=10):
#     pi2 = RS.uniform(-5, 5, n)
#     pi1 = (2/9)*(1-pi2)

#     data = np.column_stack([pi2, pi1])
#     expr = "π1=(2/9)(1-π2)"
#     save_data("stokes_pi_groups", data, expr)

# def poiseuille_pi_groups_data(n=10):
#     pi2 = RS.uniform(0, 5, n)
#     pi1 = (128.0/pi) * pi2

#     data = np.column_stack([pi2, pi1])
#     expr = "π1  = (128/π)*π2"
#     save_data("poiseuille_pi_groups", data, expr)

# def poiseuille_pressure_drop_data(n=10000):
#     # realistic fluid viscosities (water → oils)
#     mu    = RS.uniform(1e-3, 1.0,   n)      # Pa·s
#     L     = RS.uniform(0.1,  10.0,  n)      # m
#     D     = 10**RS.uniform(np.log10(0.01), np.log10(0.1), n)  # m
#     Q     = RS.uniform(1e-6,  1e-3,  n)      # m³/s
#     CONST = 128/np.pi

#     dP    = CONST * mu * L * Q * D**(-4)

#     data  = np.column_stack([mu, L, D, Q, dP])
#     # filter NaN / Inf
#     mask  = np.isfinite(data).all(axis=1)
#     dropped = n - mask.sum()
#     if dropped:
#         print(f"  [poiseuille_dp] dropped {dropped}/{n} invalid samples")
#     data  = data[mask]

#     expr = "ΔP = 128*mu*L/(π*D**4) * Q"
#     save_data("poiseuille_dp", data, expr)

# def stokes_terminal_velocity_data(n=10000):
#     # sphere densities from typical minerals → metals
#     rho_s = RS.uniform(2000, 8000,    n)    # kg/m³
#     rho_f = RS.uniform(800,  1200,    n)    # kg/m³
#     g     = RS.uniform(9.7,  9.9,     n)    # m/s²
#     r     = 10**RS.uniform(-6,  -3,   n)    # m (1 µm → 1 mm)
#     mu    = 10**RS.uniform(-5,  -2,   n)    # Pa·s

#     v_t   = (2/9) * (rho_s - rho_f) * g * r**2 / mu

#     data  = np.column_stack([rho_s, rho_f, g, r, mu, v_t])
#     mask  = np.isfinite(data).all(axis=1)
#     dropped = n - mask.sum()
#     if dropped:
#         print(f"  [stokes_vt] dropped {dropped}/{n} invalid samples")
#     data  = data[mask]

#     expr = "v_t = (2/9)*(rho_s - rho_f)*g*r**2/mu"
#     save_data("stokes_terminal_velocity", data, expr)


# def fermi_dos_with_variable_constants(n=10000):
#     V    = RS.uniform(0.1,  5.0,    n)  
#     Eev  = RS.uniform(0.1,  5.0,    n)

#     # true constants ±20%
#     ME   = 9.10938356e-31
#     HBAR = 1.054571817e-34
#     EV2J = 1.602176634e-19

#     me   = RS.uniform(0.8*ME,   1.2*ME,   n)
#     hbar = RS.uniform(0.8*HBAR, 1.2*HBAR, n)
#     ev2j = RS.uniform(0.8*EV2J, 1.2*EV2J, n)

#     E    = Eev * ev2j
#     g    = (V / (2 * np.pi**2)) * (2 * me / hbar**2)**1.5 * np.sqrt(E)

#     data = np.column_stack([V, E, me, hbar, ev2j, g])
#     mask = np.isfinite(data).all(axis=1)
#     dropped = n - mask.sum()
#     if dropped:
#         print(f"  [fermi_DOS] dropped {dropped}/{n} invalid samples")
#     data = data[mask]

#     expr = "g = V/(2π²) * (2 me/ħ²)^(3/2) * E^(1/2)"
#     save_data("fermi_gas_DOS", data, expr)


# def field_energy_with_variable_constants(n=10000):
#     E    = RS.uniform(0.1, 5.0, n)   # V/m
#     B    = RS.uniform(0.1, 5.0, n)   # T

#     # ε0 ≈ 8.85e-12 F/m; μ0 ≈ 4π×1e-7 H/m
#     eps0 = RS.uniform(7e-12, 1e-11, n)
#     mu0  = RS.uniform(1e-6,  6e-6,  n)

#     u    = 0.5 * eps0 * E**2 + 0.5 * B**2 / mu0

#     data = np.column_stack([E, B, eps0, mu0, u])
#     mask = np.isfinite(data).all(axis=1)
#     dropped = n - mask.sum()
#     if dropped:
#         print(f"  [field_energy] dropped {dropped}/{n} invalid samples")
#     data = data[mask]

#     expr = "u = 0.5*ε0*E**2 + 0.5*B**2/μ0"
#     save_data("field_energy_density", data, expr)







def poiseuille_pressure_drop_data(n=100):
    # X[:,0]=mu, X[:,1]=L, X[:,2]=Q, X[:,3]=D
    X  = np.random.rand(n, 4)
    mu, L, Q, D = X.T
    dP = (128/np.pi) * mu * L * Q / (D**4)
    data = np.column_stack([X, dP[:, np.newaxis]])
    expr = "ΔP = 128*mu*L/(π*D**4) * Q"
    save_data("poiseuille_dp_100", data, expr)

def stokes_terminal_velocity_data(n=100):
    # X[:,0]=rho_s, X[:,1]=rho_f, X[:,2]=g, X[:,3]=r, X[:,4]=mu
    X     = np.random.rand(n, 5)
    rho_s, rho_f, g, r, mu = X.T
    vt    = (2/9) * (rho_s - rho_f) * g * r**2 / mu
    data  = np.column_stack([X, vt[:, np.newaxis]])
    expr  = "v_t = (2/9)*(rho_s - rho_f)*g*r**2/mu"
    save_data("stokes_terminal_velocity_100", data, expr)

def fermi_dos_with_variable_constants(n=100):
    # X[:,0]=V, X[:,1]=Eev, X[:,2]=me, X[:,3]=hbar, X[:,4]=ev2j
    X     = np.random.rand(n, 5)
    V, Eev, me, hbar, ev2j = X.T
    E     = Eev * ev2j
    g     = (V / (2 * np.pi**2)) * (2 * me / hbar**2)**1.5 * np.sqrt(E)
    data  = np.column_stack([X, g[:, np.newaxis]])
    expr  = "g = V/(2π²) * (2 me/ħ²)^(3/2) * E^(1/2)"
    save_data("fermi_gas_DOS", data, expr)

def field_energy_with_variable_constants(n=100):
    # X[:,0]=E, X[:,1]=B, X[:,2]=eps0, X[:,3]=mu0
    X     = np.random.rand(n, 4)
    E, B, eps0, mu0 = X.T
    u     = 0.5 * eps0 * E**2 + 0.5 * B**2 / mu0
    data  = np.column_stack([X, u[:, np.newaxis]])
    expr  = "u = 0.5*ε0*E**2 + 0.5*B**2/μ0"
    save_data("field_energy_density", data, expr)

def fermi_dos(n=1000, filename="fermi_gas_DOS"):
    """
    Compute DOS for a free-electron Fermi gas.
    Only V and Eev are treated as random variables;
    me, hbar, ev2j are held constant.
    Saves columns [V, Eev, g] to `filename`.
    """
    # sample only the true variables
    V, Eev = np.random.rand(2, n)
    
    # compute energy in joules
    E = Eev * EV2J
    
    # density of states
    g = (V / (2 * np.pi**2)) * (2 * ME / HBAR**2)**1.5 * np.sqrt(E)
    
    # stack only the variable columns + result
    data = np.column_stack([V, Eev, g])
    expr = "g = V/(2π²) * (2·me/ħ²)^(3/2) · E^(1/2)"
    save_data(filename, data, expr)


def field_energy(n=1000, filename="field_energy_density"):
    """
    Compute electromagnetic field energy density.
    Only E and B are treated as random variables;
    eps0, mu0 are held constant.
    Saves columns [E, B, u] to `filename`.
    """
    # sample only the true variables
    E, B = np.random.rand(2, n)
    
    # energy density
    u = 0.5 * EPS0 * E**2 + 0.5 * B**2 / MU0
    
    # stack only the variable columns + result
    data = np.column_stack([E, B, u])
    expr = "u = 0.5·ε₀·E² + 0.5·B²/μ₀"
    save_data(filename, data, expr)

# ---------------------------------------------------------------------
if __name__ == "__main__":
    os.makedirs(".", exist_ok=True)   # use cwd

    # field_energy_original()
    # field_energy_pi()

    # fermi_dos_original()
    # fermi_dos_pi()

    # fermi_dos_with_variable_constants()
    # field_energy_with_variable_constants()
    # stokes_terminal_velocity_data()
    # poiseuille_pressure_drop_data()
    # stokes_pi_groups_data()
    # poiseuille_pi_groups_data()

    # doppler_original()
    # doppler_pi()

    field_energy()
    fermi_dos()

    print("All files written.")
