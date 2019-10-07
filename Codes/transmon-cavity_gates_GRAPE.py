# Imports
import qutip as qt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import qutip.control.pulseoptim as cpo
import qutip.logging_utils as logging
import scipy.optimize as spopt
import scipy as sp

results_folder = "/transmon-cavity_results"

# Constants
h_bar = 1

cavity_levels = 6

# Basic operators
# Atom(a and a-dagger)
a = qt.tensor(qt.destroy(2), qt.qeye(cavity_levels))
ad = a.dag()
# Cavity(c, c-dagger, and sigmaZ)
c = qt.tensor(qt.qeye(2), qt.destroy(cavity_levels))
cd = c.dag()

# Simulation time and accuracy details
num_time_steps = 5  # Number of time states for the simulation
total_time = 5  # Total simulation time
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
H_ac = g * (cd * a + c * ad) # Atom-Cavity interaction hamiltonian
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


def get_drive_fields_grape(atom_operator, cavity_operator, disp_graphs=False):
    # Target state
    targ = qt.tensor(atom_operator, cavity_operator)

    result = qt.control.grape_unitary(targ, H0, H_d, num_grape_iter, times)

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

    return QI_a, QQ_a, QI_c, QQ_c


def run_operator(state, QI_a, QQ_a, QI_c, QQ_c):
    atom_init, cavity_init = qt.ptrace(state, 0), qt.ptrace(state, 1)

    prod_atom = 1
    prod_cavity = 1
    for k in range(num_time_steps):
        H = H0 + Ha_I*QI_a[k] + Ha_Q*QQ_a[k] + Hc_I*QI_c[k] + Hc_Q*QQ_c[k]
        U_atom = ((-1j * dt) * H).expm()
        U_cavity = ((-1j * dt) * H).expm()
        prod_atom = U_atom * prod_atom
        prod_cavity = U_cavity * prod_cavity
    atom_final = qt.ptrace(prod_atom*state, 0)
    cavity_final = qt.ptrace(prod_cavity*state, 1)

    return atom_final, cavity_final


QI_a, QQ_a, QI_c, QQ_c = get_drive_fields_grape(qt.hadamard_transform(), qt.qeye(cavity_levels), disp_graphs=True)
state_init = qt.tensor(qt.basis(2,0), qt.basis(cavity_levels, 0))
#atom, cavity = run_operator(qt.tensor(qt.basis(2,0), qt.basis(cavity_levels, 0)), QI_a, QQ_a, QI_c, QQ_c)
A = [1]*len(QI_a)
state_final = qt.mesolve([H0, [Ha_I, QI_a], [Ha_Q, QQ_a], [Hc_Q, QQ_c]], [Hc_Q, QQ_c], state_init, times)
atom = qt.ptrace(state_final, 0)
bl = qt.Bloch()
bl.add_states(atom)
# bl.add_states(cavity)
bl.show()
# qt.plot_wigner_fock_distribution(qt.basis(2, 1))
