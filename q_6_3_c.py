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

'''
#Convergence delta t 
xfin_list = []
dt_list = []
dt_list_square = []

# Liste des valeurs de nsteps à tester
#nsteps_values = np.linspace(600, 2000, 10, dtype=int)
nsteps_values = [870,900,930,950,1000,1020,1050,1070,1250]

for nsteps in nsteps_values:
    dt = tfin / nsteps
    dt_list.append(dt)
    dt_list_square.append(dt*dt)
    

    
    output_name = f"conv_nsteps_{nsteps}"

    params = fct.read_in_file(input_filename)
    params.update({
        "tfin" : tfin,
        "xL" : xL,
        "xR" : xR,
        "xa" : xa,
        "xb" : xb,
        "om0" : om0,
        "V0" : V0,
        "x0" : x0,
        "n" : n,
        "sigma_norm" : sigma_norm,
        "Nsteps": nsteps,
        "Nintervals": Nintervals,
        "output": output_name
    })

    fct.run_simulation(executable, input_filename, output_name, **params)

    # Lire les résultats
    x, t, _, _, _, _, _, _, xmoy, _, _, _, _ = fct.read_quantum_data(f"outputs/{output_name}")

    # Récupérer la position moyenne à t_fin
    xfin_list.append(xmoy[-1])

plt.figure()
plt.plot(dt_list_square, xfin_list, marker='o', linestyle='-', color='blue', label=r"$\langle x \rangle(t_{fin})$")
plt.xlabel(r"$(\Delta t)^2$ ", fontsize=16)
plt.ylabel(r"$\langle x \rangle(t_{fin})$ ", fontsize=16)
plt.grid(True)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=14)
plt.tight_layout()
fct.save_figure("convergence_dt2_xfin.png")
plt.show()
'''


#Convergence h_x
xfin_list = []
dx_list = []
dx_list_square = []

# Liste des valeurs de nsteps à tester
Ninter_values = np.linspace(300, 1200, 10, dtype=int)

for ninterval in Ninter_values:
    dx = (xR - xL) / ninterval
    dx_list.append(dx)
    dx_list_square.append(dx*dx)
    
    output_name = f"conv_nx_{ninterval}"

    params = fct.read_in_file(input_filename)
    params.update({
        "tfin" : tfin,
        "xL" : xL,
        "xR" : xR,
        "xa" : xa,
        "xb" : xb,
        "om0" : om0,
        "V0" : V0,
        "x0" : x0,
        "n" : n,
        "sigma_norm" : sigma_norm,
        "Nsteps": Nsteps,
        "Nintervals": ninterval,
        "output": output_name
    })

    fct.run_simulation(executable, input_filename, output_name, **params)

    # Lire les résultats
    x, t, _, _, _, _, _, _, xmoy, _, _, _, _ = fct.read_quantum_data(f"outputs/{output_name}")

    # Récupérer la position moyenne à t_fin
    xfin_list.append(xmoy[-1])

plt.figure()
plt.plot(dx_list_square, xfin_list, marker='o', linestyle='-', color='blue', label=r"$\langle x \rangle(t_{fin})$")
plt.xlabel(r"$(\Delta x)^2$ ", fontsize=16)
plt.ylabel(r"$\langle x \rangle(t_{fin})$ ", fontsize=16)
plt.grid(True)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=14)
plt.tight_layout()
fct.save_figure("convergence_dx2_xfin.png")
plt.show()
