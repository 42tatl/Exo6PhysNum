import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 14
})

output_name = "outputs/probe_E0_psi2.out"
input_filename = "config_6_4.in"
params = fct.read_in_file(input_filename)

(tfin, xL, xR, xa, xb, om0, _, x0, n, sigma_norm, Nsteps, Nintervals) = fct.get_simulation_params(params)

dx = (xR - xL) / Nintervals

# === Load and reshape the wavefunction data ===
raw_data = np.loadtxt(output_name)

# Reshape to (Nsteps, Npoints * 3) for abs, real, imag
if raw_data.ndim == 1:
    raw_data = raw_data[np.newaxis, :]  # ensure at least 2D

Nsteps = raw_data.shape[0]

# Extract Re and Im
re_psi = raw_data[:, 1::3]
im_psi = raw_data[:, 2::3]

# Compute norm at each time step: ∫ |ψ|² dx = ∑ (Re² + Im²) * dx
norm_t = np.sum(re_psi**2 + im_psi**2, axis=1) * dx

# === Plot norm over time ===
plt.figure(figsize=(7, 4))
plt.plot(norm_t, label=r"$\|\psi(t)\|^2$")
plt.axhline(1.0, color='red', linestyle='--', label=r"Expected = 1")
plt.xlabel(r"$nsteps$")
plt.ylabel(r"$\int |\psi|^2 dx$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/norm_conservation_check.pdf")
plt.show()
