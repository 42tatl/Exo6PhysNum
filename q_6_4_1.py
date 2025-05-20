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

executable = './Exe.exe'
repertoire = r"C:\\Users\\Avril\\Desktop\\Exo6PhysNum"
os.chdir(repertoire)

input_filename = "config_6_4.in"
params = fct.read_in_file(input_filename)

(tfin, xL, xR, xa, xb, 
om0, V0, x0, n, sigma_norm, 
Nsteps, Nintervals, output, verbose) = fct.get_simulation_params(params)

os.makedirs("outputs", exist_ok=True)

# === Réponse à la question (i) ===
t_trans = 0.035
V0_values = [0.1, 0.5, 1.0, 2.0, 5.0]

E_over_V0 = []
P_right_list = []

for V0 in V0_values:
    output_name = f"transmission_V0_{V0}"
    params["V0"] = V0
    params["output"] = output_name

    success = fct.run_simulation(executable, input_filename, output_name, **params)
    if not success:
        print(f"Simulation failed for V0 = {V0}")
        continue

    x, t, psi2, obs = fct.read_quantum_data(f"outputs/{output_name}")

    idx_trans = np.argmin(np.abs(t - t_trans))
    E_t = obs[idx_trans, 3]         # Energie moyenne
    P_right = obs[idx_trans, 2]     # Probabilité x > 0

    E_over_V0.append(E_t / V0)
    P_right_list.append(P_right)

# Tracer la courbe demandée
plt.figure(figsize=(8, 5))
plt.plot(E_over_V0, P_right_list, 'o-', label=r"$P_{x>0}(t_{\mathrm{trans}})$")
plt.xlabel(r"$\langle E \rangle / V_0$")
plt.ylabel(r"Probabilit\'e de transmission $P_{x>0}$")
plt.title("Transmission en fonction de l'\u00e9nergie relative")
plt.grid(True)
plt.legend()
fct.save_figure("transmission_vs_EsurV0.pdf")
plt.show()