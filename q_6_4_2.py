import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct
from matplotlib.animation import FuncAnimation

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 13
})

executable = './Exe.exe'
repertoire = r"C:\\Users\\Avril\\Desktop\\Exo6PhysNum"
os.chdir(repertoire)

input_filename = "config_6_4.in"
params = fct.read_in_file(input_filename)

(tfin, xL, xR, xa, xb, 
om0, V0, x0, n, sigma_norm, 
Nsteps, Nintervals) = fct.get_simulation_params(params)

os.makedirs("outputs", exist_ok=True)

# === Réponse à la question (ii) ===
# On choisit 3 cas : V0 < E, V0 ~ E, V0 > E
V0_test = [0.1, 0.4, 1.0]  # À ajuster après test pour approximer V0/E
labels = ["V0 < E", "V0 ~ E", "V0 > E"]
titles = [r"$V_0 < \langle E \rangle$", r"$V_0 \approx \langle E \rangle$", r"$V_0 > \langle E \rangle$"]

for V0, label, title in zip(V0_test, labels, titles):
    output_name = f"evolution_{label.replace(' ', '_')}"
    params["V0"] = V0
    params["output"] = output_name

    success = fct.run_simulation(executable, input_filename, output_name, **params)
    if not success:
        print(f"Simulation failed for V0 = {V0}")
        continue

    x, t, psi2, obs = fct.read_quantum_data(f"outputs/{output_name}")

    # === Tracer |psi(x,t)|^2 à différents temps ===
    times_to_plot = [0, int(len(t) * 0.25), int(len(t) * 0.5), int(len(t) * 0.75), -1]
    plt.figure(figsize=(8, 5))
    for idx in times_to_plot:
        plt.plot(x, psi2[idx], label=fr"$t = {t[idx]:.3f}$")
    plt.xlabel("Position x")
    plt.ylabel(r"$|\psi(x,t)|^2$")
    plt.title(title)
    plt.grid()
    plt.legend()
    fct.save_figure(f"psi2_snapshots_{label.replace(' ', '_')}.pdf")
    plt.show()

    # === Tracer P_{x<0}(t) et P_{x>0}(t) ===
    Px_left = obs[:, 1]
    Px_right = obs[:, 2]

    plt.figure(figsize=(7, 4))
    plt.plot(t, Px_left, label=r"$P_{x<0}(t)$")
    plt.plot(t, Px_right, label=r"$P_{x>0}(t)$")
    plt.xlabel("Temps t")
    plt.ylabel("Probabilité")
    plt.title(f"{title} : évolution des probabilités")
    plt.grid(True)
    plt.legend()
    fct.save_figure(f"probas_{label.replace(' ', '_')}.pdf")
    plt.show()