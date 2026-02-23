#!/usr/bin/env python3
# Head-on binary droplet coalescence — We = 2, Re ~ 236.
#
# Two water droplets (D = 306 um) in air collide along the x-axis.
# Ported from an OpenFOAM interFoam simulation (temp/ directory).
#
# Nondimensional parameters follow Meng (2016) convention:
#   rho_l = 100, rho_g = 1, sigma = 0.72
#   gamma_l = 4.4, pi_inf_l = 100;  gamma_g = 1.4, pi_inf_g = 0
#   D = 1, R = 0.5, p_gas = 1.0, p_liq = 1 + sigma/R = 2.44
#   Ur = sqrt(We * sigma / (rho_l * D)) = 0.12
#   Re = rho_l * Ur * D / mu_l = 236
#
# 3 patches: gas background + two spherical droplets (no hardcoded IC needed).

import math
import json

# --- Nondimensional parameters (Meng 2016) ---
rho_l = 100.0
rho_g = 1.0
sigma = 0.72
gamma_l, pi_inf_l = 4.4, 100.0
gamma_g, pi_inf_g = 1.4, 0.0
p_gas = 1.0

D = 1.0
R = D / 2.0
p_liq = p_gas + sigma / R  # Laplace pressure jump

# Collision parameters
We = 2.0
Re = 236.0
Ur = math.sqrt(We * sigma / (rho_l * D))  # ~0.12 relative velocity
Ud = Ur / 2.0  # each droplet moves at half the relative velocity

# Viscosity (matching Re = 236)
mu_l = rho_l * Ur * D / Re   # ~0.05085
mu_g = mu_l / 48.1           # ~0.001057 (physical viscosity ratio water/air)

# Domain: 9D x 6D x 6D centered at origin
x0, x1 = -4.5, 4.5
y0, y1 = -3.0, 3.0
z0, z1 = -3.0, 3.0

# Grid: 10 cells/D
Nx = 179   # 90 cells in x (9D * 10)
Ny = 119   # 60 cells in y (6D * 10)
Nz = 119   # 60 cells in z (6D * 10)

eps = 1e-9

# Droplet placement: centers at (+-sep, 0, 0)
# Minimal gap (~0.3D) to avoid wake-induced deformation before contact
sep = 0.65 * D  # center-to-center = 1.3D, surface gap = 0.3D

# Time stepping
dx = (x1 - x0) / (Nx + 1)
c_l = math.sqrt(gamma_l * (p_gas + pi_inf_l) / rho_l)  # ~2.11
mydt = 0.1 * dx / c_l  # CFL = 0.1

# Collision timescale = D / Ur ~ 8.33
# Run for ~11 collision times (matching 5ms OpenFOAM end time)
t_end = 11.0 * D / Ur  # ~92
Nt = int(math.ceil(t_end / mydt))
Ns = max(1, Nt // 100)  # ~100 output frames

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "z_domain%beg": z0,
    "z_domain%end": z1,
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Ns,
    # Simulation Algorithm
    "model_eqns": 3,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "F",
    "time_stepper": 3,
    "weno_order": 3,
    "avg_state": 2,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "T",
    "bc_x%beg": -6,
    "bc_x%end": -6,
    "bc_y%beg": -6,
    "bc_y%end": -6,
    "bc_z%beg": -6,
    "bc_z%end": -6,
    "num_patches": 3,
    "num_fluids": 2,
    "weno_avg": "T",
    "surface_tension": "T",
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "alpha_wrt(1)": "T",
    "cf_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "vel_wrt(3)": "T",
    "parallel_io": "T",
    "sigma": sigma,
    # Fluid Parameters (Liquid — fluid 1)
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    # Fluid Parameters (Gas — fluid 2)
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu_g,
    # ===== Patch 1: Gas background (fills entire domain) =====
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": x1 - x0,
    "patch_icpp(1)%length_y": y1 - y0,
    "patch_icpp(1)%length_z": z1 - z0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val": 0,
    # ===== Patch 2: Left droplet (center at -sep, moves in +x) =====
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": -sep,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.95,
    "patch_icpp(2)%vel(1)": Ud,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p_liq,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1,
    # ===== Patch 3: Right droplet (center at +sep, moves in -x) =====
    "patch_icpp(3)%geometry": 8,
    "patch_icpp(3)%x_centroid": sep,
    "patch_icpp(3)%y_centroid": 0.0,
    "patch_icpp(3)%z_centroid": 0.0,
    "patch_icpp(3)%radius": R,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%smoothen": "T",
    "patch_icpp(3)%smooth_patch_id": 1,
    "patch_icpp(3)%smooth_coeff": 0.95,
    "patch_icpp(3)%vel(1)": -Ud,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%vel(3)": 0.0,
    "patch_icpp(3)%pres": p_liq,
    "patch_icpp(3)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(3)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(3)%alpha(1)": 1.0 - eps,
    "patch_icpp(3)%alpha(2)": eps,
    "patch_icpp(3)%cf_val": 1,
}

print(json.dumps(data))
