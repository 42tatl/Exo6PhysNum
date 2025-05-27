import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct


executable = './Exe'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo6PhysNum"  # Modify for correct directory
os.chdir(repertoire)
output_name = "transmission_6_3"

input_filename = "config_6_3.in"
params = fct.read_in_file(input_filename)

(tfin, xL, xR, xa, xb, 
om0, V0, x0, n, sigma_norm, Nsteps, Nintervals) = fct.get_simulation_params(params)

fct.run_simulation(executable, input_filename, output_name, **params)
x, t, abs_psi, re_psi, im_psi, P_left, P_right, E, xmoy, x2moy, pmoy, p2moy, uncertainty = fct.read_quantum_data(f"outputs/{output_name}")
hbar = 1

P_total = P_left + P_right

E_mean = np.mean(E)
delta_E = np.abs(np.diff(E))  # longueur N-1
t_mid = t[:-1] 

plt.figure()
plt.plot(t, P_total, label="Total probability")
plt.xlabel(r"t", fontsize=16)
plt.ylabel(r"Total Probability", fontsize=16)
plt.ylim(0.98, 1.02)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.grid()
plt.legend(fontsize=14)
fct.save_figure("proba_totale.png")
plt.show()

plt.figure()
plt.plot(t, E, label=r"⟨H(t)⟩")
plt.axhline(y=E_mean, color='red', linestyle='--', label=fr"⟨H(0)⟩ = {E_mean:.3e}")
plt.xlabel(r"t ", fontsize=16)
plt.ylabel(r"Energy ", fontsize=16)
plt.ylim(2000,3000)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.grid()
plt.legend(fontsize=14)
fct.save_figure("energie_moyenne_vs_t.png")
plt.show()


plt.figure()
plt.plot(t, uncertainty, label=r"$⟨\Delta x(t)⟩ \cdot ⟨\Delta p(t)⟩$")
plt.axhline(y=hbar/2, color='red', linestyle='--', label=r"$\hbar/2$")
plt.xlabel(r"t ", fontsize=16)
plt.ylabel(r"$⟨\Delta x(t)⟩ \cdot ⟨\Delta p(t)⟩$", fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.grid()
plt.legend(fontsize=14)
fct.save_figure("incertitude_vs_t.png")
plt.show()