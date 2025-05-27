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

executable = './Exe'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo6PhysNum"  # Modify for correct directory
os.chdir(repertoire)
output_name = "transmission_6_3"

input_filename = "config_6_3.in"
params = fct.read_in_file(input_filename)

(tfin, xL, xR, xa, xb, 
om0, V0, x0, n, sigma_norm, Nsteps, Nintervals) = fct.get_simulation_params(params)

fct.run_simulation(executable, input_filename, output_name, **params)
x, t, abs_psi, re_psi, im_psi, P_left, P_right, E, xmoy, x2moy, pmoy, p2moy = fct.read_quantum_data(f"outputs/{output_name}")

#QUESTION A

#Analytique
print(pmoy[0])
print(xmoy.shape)
print(pmoy.shape)
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

if x is not None:
    f = np.array(abs_psi)  # f = |ψ(x,t)|

    # On veut f[x, t] (colonne = x, ligne = t), donc on transpose si nécessaire
    f_plot = f.T if f.shape[0] == len(t) else f

    
    im = plt.imshow(f_plot, aspect='auto',
                    extent=[x[0], x[-1], t[0], t[-1]],
                    origin='lower', cmap='turbo')

    cbar = plt.colorbar(im)
    cbar.set_label(r"$|\psi(x,t)|$", fontsize=16)

    plt.ylabel(r"$t$ [s]", fontsize=16)
    plt.xlabel(r"$x$ [m]", fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    fct.save_figure("psi_abs_t_vs_x.png")
    plt.show()

if x is not None:
    f = np.array(re_psi)  # f = |ψ(x,t)|

    # On veut f[x, t] (colonne = x, ligne = t), donc on transpose si nécessaire
    f_plot = f.T if f.shape[0] == len(t) else f

    
    im = plt.imshow(f_plot, aspect='auto',
                    extent=[x[0], x[-1], t[0], t[-1]],
                    origin='lower', cmap='turbo')

    cbar = plt.colorbar(im)
    cbar.set_label(r"$Re(\psi(x,t))$", fontsize=16)

    plt.ylabel(r"$t$ [s]", fontsize=16)
    plt.xlabel(r"$x$ [m]", fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    fct.save_figure("psi_Re_t_vs_x.png")
    plt.show()

