import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from matplotlib.animation import FuncAnimation

def read_in_file(filename):
    variables = {}
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("#") and '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.split('!')[0].strip()
                if value.lower() in ('true', 'false'):
                    variables[key] = value.lower() == 'true'
                elif value.replace('.', '').replace('e-', '').replace('e+', '').isdigit():
                    if '.' in value or 'e' in value.lower():
                        variables[key] = float(value)
                    else:
                        variables[key] = int(value)
                else:
                    variables[key] = value
    return variables

def run_simulation(executable, input_filename, output_name, **params):
    '''Runs the simulation with the given parameters'''
    os.makedirs("outputs", exist_ok=True)

    # Build command
    cmd = [executable, input_filename]
    for key, value in params.items():
        if key != "output":
            cmd.append(f"{key}={value}")
    cmd.append(f"output=outputs/{output_name}")

    print("\nRunning command:", " ".join(cmd))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        print("Command executed successfully")

        # Check for expected output files (.out extensions used in quantum case)
        expected_files = [
            f"outputs/{output_name}_pot.out",
            f"outputs/{output_name}_psi2.out",
            f"outputs/{output_name}_obs.out"
        ]

        if all(os.path.exists(f) for f in expected_files):
            print("All output files created successfully")
            return True
        else:
            missing = [f for f in expected_files if not os.path.exists(f)]
            print(f"Missing files: {missing}")
            return False

    except subprocess.CalledProcessError as e:
        print(f"Command failed with error:\n{e.stderr}")
        return False
    
def get_simulation_params(params):
    '''Extracts all quantum simulation parameters with proper typing'''

    # Physical parameters
    tfin = float(params.get("tfin", 0.08))
    xL = float(params.get("xL", -1.0))
    xR = float(params.get("xR", 1.0))
    xa = float(params.get("xa", -0.5))
    xb = float(params.get("xb", 0.5))
    om0 = float(params.get("om0", 100.0))
    V0 = float(params.get("V0", 0.0))
    x0 = float(params.get("x0", -0.5))
    n = int(params.get("n", 16))
    sigma_norm = float(params.get("sigma_norm", 0.04))

    # Numerical parameters
    Nsteps = int(params.get("Nsteps", 800))
    Nintervals = int(params.get("Nintervals", 512))

    return (
        tfin, xL, xR, xa, xb, om0, V0,
        x0, n, sigma_norm, Nsteps, Nintervals
    )

def read_quantum_data(output_base):
    """
    Lit les fichiers de sortie de la simulation quantique.
    Retourne : x, t, psi2, obs
    """
    try:
        # Lecture de la grille x depuis _pot.out (colonne 0)
        pot_data = np.loadtxt(f"{output_base}_pot.out")
        if pot_data.ndim == 1:
            x = pot_data[0:1]  # cas bord (1 ligne)
        else:
            x = pot_data[:, 0]

        # Lecture de psi^2 (chaque ligne = 1 pas de temps)
        psi2 = np.loadtxt(f"{output_base}_psi2.out")
        if psi2.ndim == 1:
            psi2 = psi2.reshape((1, -1))

        # Lecture des observables (col. 0 = t)
        obs = np.loadtxt(f"{output_base}_obs.out")
        if obs.ndim == 1:
            obs = obs.reshape((1, -1))

        t = obs[:, 0]
        P_left = obs[:, 1]   # Probabilité x < 0
        P_right = obs[:, 2]  # Probabilité x > 0
        E = obs[:, 3]        # Energie moyenne
        xmoy = obs[:, 4]     # Position moyenne
        x2moy = obs[:, 5]    # Position carrée moyenne
        pmoy = obs[:, 6]     # Quantité de mouvement moyenne
        p2moy = obs[:, 7]    # Quantité de mouvement carrée moyenne

        psi_raw = np.loadtxt(f"{output_base}_psi2.out")
        Npoints = len(x)
        Nsteps = len(t)

        psi_raw = psi_raw.reshape((Nsteps, Npoints * 3))

        abs_psi = psi_raw[:, 0::3]
        re_psi = psi_raw[:, 1::3]  # => Re(ψ)
        im_psi = psi_raw[:, 2::3]  # => Im(ψ)
        return x, t, abs_psi, re_psi, im_psi, P_left, P_right, E, xmoy, x2moy, pmoy, p2moy
    except Exception as e:
        print(f"Error reading quantum data: {e}")
        return None, None, None, None, None, None, None, None, None, None
    




def animate_quantum_psi(x, t, psi2, save_as=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    line, = ax.plot(x, psi2[0, :], 'b-')
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, np.max(psi2))
    ax.set_xlabel('Position x')
    ax.set_ylabel(r'$|\psi(x,t)|^2$')
    ax.set_title('Quantum Probability Density Evolution')
    ax.grid(True)

    def update(frame):
        line.set_ydata(psi2[frame, :])
        ax.set_title(f"t = {t[frame]:.3f}")
        return line,

    ani = FuncAnimation(fig, update, frames=len(t), interval=50, blit=True)

    if save_as:
        try:
            ani.save(save_as, writer='ffmpeg', fps=30)
        except:
            ani.save(save_as.replace('.mp4', '.gif'), writer='pillow', fps=15)

    plt.show()
    return ani

def save_figure(filename, fig=None, subfolder="figures", dpi=300, tight=True):
    if fig is None:
        fig = plt.gcf()
    os.makedirs(subfolder, exist_ok=True)
    filepath = os.path.join(subfolder, filename)
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight' if tight else None)
    print(f"Figure saved to {filepath}")

def plot_observables(obs):
    t = obs[:, 0]
    E = obs[:, 3]
    xm = obs[:, 4]
    pm = obs[:, 6]

    plt.figure()
    plt.plot(t, E, label="Energie")
    plt.plot(t, xm, label=r"$\langle x \rangle$")
    plt.plot(t, pm, label=r"$\langle p \rangle$")
    plt.xlabel("Temps")
    plt.legend()
    plt.title("\u00c9volution des observables")
    plt.grid()
    save_figure("observables.pdf")
    plt.show()

# Example usage:
# output_base = "outputs/simulation1"
# x, t, psi2, obs = read_quantum_data(output_base)
# animate_quantum_psi(x, t, psi2, save_as="quantum_evolution.mp4")
# plot_observables(obs)