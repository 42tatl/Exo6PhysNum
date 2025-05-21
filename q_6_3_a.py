import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct

#QUESTION A
executable = './Exe2'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo6PhysNum"  # Modify for correct directory
os.chdir(repertoire)
output_name = "transmission_6_3"

input_filename = "config_6_3.in"
params = fct.read_in_file(input_filename)

(tfin, xL, xR, xa, xb, 
om0, V0, x0, n, sigma_norm, Nsteps, Nintervals) = fct.get_simulation_params(params)

fct.run_simulation(executable, input_filename, output_name, **params)
x, t, psi2, P_left, P_right, E, xmoy, x2moy, pmoy, p2moy = fct.read_quantum_data(f"outputs/{output_name}")


#Analytique
print(pmoy[0])
xclass = x0 * np.cos(om0 * t) + pmoy[0] * np.sin(om0 * t) / om0
pclass = pmoy[0] * np.cos(om0 * t) - x0 * om0 * np.sin(om0 * t)


plt.figure()
plt.plot(t, xmoy, label=r"$\langle x(t) \rangle$")
plt.plot(t, xclass, '--', label=r"$x_{\mathrm{class}}(t)$")
plt.xlabel(r"t [s]")
plt.ylabel(r"<x(t)> [m]")
plt.grid()
plt.legend()
fct.save_figure("xmoy_t.png")
plt.show()

# Tracer <p>(t)
plt.figure()
plt.plot(t, pmoy, label=r"$\langle p(t) \rangle$")
plt.plot(t, pclass, '--', label=r"$p_{\mathrm{class}}(t)$")
plt.xlabel(r"t [s]")
plt.ylabel(r"<p(t)> [kg.m/s]")
plt.grid()
plt.legend()
fct.save_figure("pmoy_t.png")
plt.show()



