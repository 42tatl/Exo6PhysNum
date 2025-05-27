import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize


# Plot style
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 14
})

# Paths
repertoire = r"C:/Users/Avril/Desktop/Exo6PhysNum"
executable = os.path.normpath(os.path.join(repertoire, "Exe.exe"))
os.chdir(repertoire)

# Read input
input_filename = "config_6_4.in"
params = fct.read_in_file(input_filename)
(tfin, xL, xR, xa, xb, om0, V0, x0, n, sigma_norm, Nsteps, Nintervals) = fct.get_simulation_params(params)
os.makedirs("outputs", exist_ok=True)

# --- Estimate initial mean energy ----------------------------------
probe_name = "probe_E0"
params["V0"] = 0.0
params["output"] = probe_name
fct.run_simulation(executable, input_filename, probe_name, **params)
_, t, _, _, _, _, _, E_probe, _, _, _, _ = fct.read_quantum_data(f"outputs/{probe_name}")
E_mean = E_probe[0]
print(f"Estimated mean energy ⟨E⟩ = {E_mean:.3f}")

if E_probe is not None:
    plt.figure(figsize=(7, 4))
    plt.plot(t, E_probe - E_probe[0], label=r"$\Delta E$")
    plt.axhline(0, linestyle='--', color='red', label = "Expected = 0")
    plt.ylabel(r"$\langle E(t) \rangle - \langle E(0) \rangle$")
    plt.xlabel(r"$t$ [s]")
    plt.legend()
    plt.grid(True)
    fct.save_figure("64_energy_conservation_V0_0.pdf")
    plt.show()
    
else:
    print("Failed to load energy data for V0 = 0.")

# Barrier test cases relative to ⟨E⟩
V0_test = [0.5 * E_mean, 1.0 * E_mean, 2.0 * E_mean]
print(E_probe[0], V0_test)
plot_titles = [
    r"$V_0 < \langle E \rangle$",
    r"$V_0 \approx \langle E \rangle$",
    r"$V_0 > \langle E \rangle$"
]
safe_labels = ["V0_lt_E", "V0_eq_E", "V0_gt_E"]

# === Loop through the 3 barrier cases ===
for V0, label_plot, safe_label in zip(V0_test, plot_titles, safe_labels):
    output_name = f"evolution_{safe_label}"
    params["V0"] = V0
    params["output"] = output_name

    success = fct.run_simulation(executable, input_filename, output_name, **params)
    if not success:
        print(f"Simulation failed for V0 = {V0}")
        continue

    # Load wavefunction and observable data
    x, t, abs_psi, re_psi, im_psi, P_left, P_right, E, xmoy, x2moy, pmoy, p2moy = fct.read_quantum_data(f"outputs/{output_name}")
    if t is None or abs_psi is None:
        print(f"Data loading failed for V0 = {V0}")
        continue

    # === Plot 2D contour of |\psi(x,t)| ===
    
    X, T = np.meshgrid(x, t)
    plt.figure(figsize=(5,4))

    contour = plt.contourf(X, T, abs_psi, levels=100, cmap='turbo', norm=Normalize(vmin=0, vmax=3.0))

    plt.xlabel(r"$x$ [m]")  # Use [a.u.] if units are arbitrary
    plt.ylabel(r"$t$ [s]")
    plt.colorbar(contour, label=r"$|\psi(x,t)|$ [a.u.]")

    plt.tight_layout()
    fct.save_figure(f"64_psi_{safe_label}.pdf")
    plt.show()

    # === Plot P_{x<0}(t) and P_{x>0}(t) ===
    plt.figure(figsize=(5, 4))
    plt.plot(t, P_left, label=r"$P_{x<0}(t)$")
    plt.plot(t, P_right, label=r"$P_{x>0}(t)$")
    plt.xlabel(r"$t$ [s]")
    plt.ylabel(r"Probability $P(t)$")
    plt.grid(True)
    plt.legend()
    fct.save_figure(f"64_probas_{safe_label}.pdf")
    plt.show()
