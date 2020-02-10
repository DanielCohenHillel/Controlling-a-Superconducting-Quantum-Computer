# Imports
import qutip as qt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import qutip.control.pulseoptim as cpo
import qutip.logging_utils as logging
import scipy.optimize as spopt
import scipy as sp
import time

print("-------------------------------------------------------------------------------------------------------------\n")

results_folder = "/transmon-cavity_results"

# Constants
h_bar = 1

cavity_levels = 8

# Basic operators
# Atom(a and a-dagger)
a = qt.tensor(qt.destroy(2), qt.qeye(cavity_levels))
ad = a.dag()
# Cavity(c, c-dagger, and sigmaZ)
c = qt.tensor(qt.qeye(2), qt.destroy(cavity_levels))
cd = c.dag()

# Simulation time and accuracy details
num_time_steps = 10  # Number of time states for the simulation
total_time = 20  # Total simulation time
dt = total_time/num_time_steps  # Time step size
times = np.linspace(0.0, total_time, num_time_steps)  # Time space array
num_grape_iter = 10
# Frequencies
w_a = 1  # Atom(qubit) frequency
w_c = 1  # Cavity frequency
g = 1  # Linear photon number dependent dispersive shifts

# Hamiltonians (Linear/Harmonic)
H_a = w_a*a*ad  # Atom hamiltonian
H_c = w_c*c*cd  # Cavity hamiltonian
H_ac = g * (cd * c * ad * a)  # Atom-Cavity interaction hamiltonian
H_ac = g * ad*a*cd*c  # Atom-Cavity interaction hamiltonian
H0 = H_a + H_c + H_ac  # Full atom-cavity hamiltonian w/o a derive
# Drive hamiltonians
# Atom
Ha_I = a + ad
Ha_Q = 1j*(a-ad)
# Cavity
Hc_I = c + cd
Hc_Q = 1j*(c-cd)
# Time dependent(controllable) derive hamiltonian list
H_d = [Ha_I, Ha_Q, Hc_I, Hc_Q]

bl = qt.Bloch()


def get_drive_fields_grape(operator, disp_graphs=False):
    start_time = time.time_ns()

    U_start = np.ones([4, num_time_steps])/7.0
    result = qt.control.grape_unitary(
        operator, H0, H_d, num_grape_iter, times, u_start=U_start)

    # Atom(transmon)
    QI_a = result.u[-1, 0, :]
    QQ_a = result.u[-1, 1, :]
    # Cavity
    QI_c = result.u[-1, 2, :]
    QQ_c = result.u[-1, 3, :]

    # Display drive graphs
    if disp_graphs:
        fig, (ax1, ax2) = plt.subplots(2, 1)

        ax1.set_title("Transmon Drive Fields")
        ax1.plot(times, QI_a, 'r')
        ax1.plot(times, QQ_a, 'r--')
        ax1.legend(("I", "Q"))
        ax1.set_xlabel('Time (sec)')
        ax1.set_ylabel('Amplitude')

        ax2.set_title("Cavity Drive Fields")
        ax2.plot(times, QI_c, 'b')
        ax2.plot(times, QQ_c, "b--")
        ax2.legend(("I", "Q"))
        ax2.set_xlabel('Time (sec)')
        ax2.set_ylabel('Amplitude')

    print("-    Evaluating the control pulses took",
          str((time.time_ns() - start_time)*1e-9)[0:4], "Seconds")
    return QI_a, QQ_a, QI_c, QQ_c


def run_operator(state, QI_a, QQ_a, QI_c, QQ_c, show_fid_graph=False):
    start_time = time.time_ns()

    fid = []

    prod = 1
    for k in range(num_time_steps):
        H = H0 + Ha_I*QI_a[k] + Ha_Q*QQ_a[k] + Hc_I*QI_c[k] + Hc_Q*QQ_c[k]
        U = ((-1j * dt) * H).expm()
        prod = U * prod
        bl.add_states(qt.ptrace(prod * state, 0))
        fid.append(qt.fidelity(prod*state, state_target*state_init))

    result = prod * state

    if show_fid_graph:
        fig, ax = plt.subplots(1, 1)
        ax.plot(times, fid)
        ax.set_title("Fidelity Over Time")
        ax.set_xlabel('Time (sec)')
        ax.set_ylabel('Fidelity')

    print("-    Simulating the system took",
          str((time.time_ns() - start_time) * 1e-9)[0:4], "Seconds")
    return result


state_init = qt.tensor(qt.basis(2, 0), qt.basis(cavity_levels, 0))
state_target = qt.tensor(qt.qeye(2), qt.create(cavity_levels, 1))

QI_a, QQ_a, QI_c, QQ_c = get_drive_fields_grape(state_target, disp_graphs=True)

state_final = run_operator(state_init, QI_a, QQ_a,
                           QI_c, QQ_c, show_fid_graph=True)

print("-    Final fidelity:",
      str(qt.fidelity(state_target, state_final)*100)[0:9] + "%")

bl.add_states(qt.ptrace(state_init, 0))
bl.add_states(qt.ptrace(state_final, 0))
bl.show()

xvec = np.linspace(-5, 5, 500)
W_init = qt.wigner(qt.ptrace(state_init, 1), xvec, xvec)
W_final = qt.wigner(qt.ptrace(state_final, 1), xvec, xvec)
W_target = qt.wigner(qt.ptrace(state_target*state_init, 1), xvec, xvec)

fig, axes = plt.subplots(1, 3)

axes[0].contourf(xvec, xvec, W_init, 100)
axes[0].set_title("Initial")

axes[1].contourf(xvec, xvec, W_final, 100)
axes[1].set_title("Final")

axes[2].contourf(xvec, xvec, W_target, 100)
axes[2].set_title("Target")

fig.tight_layout()
plt.show()

# print("\n--\nInit: " + str(state_init) + "\n")
# print("Final: " + str(state_final) + "\n--\n")

print("\n-------------------------------------------------------------------------------------------------------------")
