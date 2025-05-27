import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct

# === Plot settings ===
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 14
})

# === Directory and executable path ===
repertoire = r"C:/Users/Avril/Desktop/Exo6PhysNum"
executable = os.path.normpath(os.path.join(repertoire, "Exe.exe"))

# Confirm the executable path is valid
print("Executable path:", executable)
print("Executable exists?", os.path.exists(executable))

# Change working directory
os.chdir(repertoire)

# === Read input and simulation parameters ===
input_filename = "config_6_4.in"
params = fct.read_in_file(input_filename)

(tfin, xL, xR, xa, xb, om0, _, x0, n, sigma_norm, Nsteps, Nintervals) = fct.get_simulation_params(params)

# === Transmission probability study ===
t_trans = 0.035
V0_values = np.linspace(325, 4000, 50)
print("V0 values:", V0_values)
E_over_V0 = []
P_right_list = []

for V0 in V0_values:
    output_name = f"6_4_V0_{int(V0)}"
    params["V0"] = V0
    params["output"] = output_name

    success = fct.run_simulation(executable, input_filename, output_name, **params)
    if not success:
        print(f"Simulation failed for V0 = {V0}")
        continue

    x, t, abs_psi, re_psi, im_psi, P_left, P_right, E, xmoy, x2moy, pmoy, p2moy = fct.read_quantum_data(f"outputs/{output_name}")

    if t is None or E is None or P_right is None:
        print(f"Data loading failed for V0 = {V0}")
        continue

    idx_trans = np.argmin(np.abs(t - t_trans))
    E_t = E[idx_trans]
    P_right_val = P_right[idx_trans]

    E_over_V0.append(E_t / V0)
    P_right_list.append(P_right_val)

# === Plot results ===
plt.figure(figsize=(8, 5))
plt.plot(E_over_V0, P_right_list, '+-', label=r"$P_{x>0}(t_{\mathrm{trans}})$")
plt.xlabel(r"$\langle E \rangle / V_0$")
plt.ylabel(r"$P_{x>0}$")
plt.grid(True)
plt.legend()
fct.save_figure("64_EsurV0.pdf")
plt.show()

