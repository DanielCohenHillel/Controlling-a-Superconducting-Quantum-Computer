# Imports
from qutip import *
from qutip.control import *
import qutip
import numpy as np
from matplotlib import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots
import qutip.control.pulseoptim as cpo
from qutip.qip import hadamard_transform
import qutip.logging_utils as logging
logger = logging.get_logger()
log_level = logging.INFO

# Basic operators
q = destroy(2)
qd = q.dag()

# Initial state of the qubit
psi0 = basis(2, 0)

# Number of time division
Ns = 20
# Time of simulation evolution
T = 1
dt = T/Ns

# Time space
times = np.linspace(0.0, T, Ns)


def main():
    H0 = setH0()
    Hq_I = q + qd
    Hq_Q = 1j*(q - qd)

    H = [H0, [Hq_I, qI], [Hq_Q, qQ]]
    psitarg = settarg()

    # Time-Independent hamiltonian H0
    H_d = H0
    # Time-Dependent hamiltonian of the drives, ~~ H = H_d + u1(t)*Hq_I + u2(t)*Hq_Q ~~
    H_c = [Hq_I, Hq_Q]
    # Fidelity error target
    fid_err_targ = 1e-10
    # Maximum iterations for the optisation algorithm
    max_iter = 200
    # Maximum (elapsed) time allowed in seconds
    max_wall_time = 120
    # Minimum gradient (sum of gradients squared)
    # as this tends to 0 -> local minima has been found
    min_grad = 1e-20
    # Type of initial guess of the drives
    p_type = 'SAW'
    # Save the results to a text file
    f_ext = "{}_n_ts{}_ptype{}.txt".format("Hi", Ns, p_type)

    # The function that does the GRAPE algorithm
    result = cpo.optimize_pulse(H_d, H_c, psi0, psitarg, Ns, T,
                                fid_err_targ=fid_err_targ, min_grad=min_grad,
                                max_iter=max_iter, max_wall_time=max_wall_time,
                                out_file_ext=f_ext, init_pulse_type=p_type,
                                log_level=log_level, gen_stats=True)

    # Plot
    fig1 = plt.figure()
    # Initial
    ax1 = fig1.add_subplot(2, 1, 1)
    ax1.set_title("Initial control amps")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Control amplitude")
    ax1.step(result.time,
             np.hstack((result.initial_amps[:, 0],
                        result.initial_amps[-1, 0])),
             where='post')
    ax1.step(result.time,
             np.hstack((result.initial_amps[:, 1],
                        result.initial_amps[-1, 0])),
             where='post')

    # Final
    ax2 = fig1.add_subplot(2, 1, 2)
    ax2.set_title("Optimised Control Sequences")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Control amplitude")
    ax2.step(result.time,
             np.hstack((result.final_amps[:, 0], result.final_amps[-1, 0])),
             where='post')
    ax2.step(result.time,
             np.hstack((result.final_amps[:, 1], result.final_amps[-1, 0])),
             where='post')

    plt.tight_layout()
    plt.show()

    # Reports the results
    result.stats.report()
    print("Final evolution\n{}\n".format(result.evo_full_final))
    print("********* Summary *****************")
    print("Final fidelity error {}".format(result.fid_err))
    print("Final gradient normal {}".format(result.grad_norm_final))
    print("Terminated due to {}".format(result.termination_reason))
    print("Number of iterations {}".format(result.num_iter))

    # Displays what the Qutip built-in grape algorithm thinks he managed to achive with the drives
    b = Bloch()
    b.add_states(result.evo_full_final)
    b.show()

    # Set results into variables
    QI = result.final_amps[:, 0]
    QQ = result.final_amps[:, 1]

    # Calls our own simulation of the qubit with the given drives from the built in GRAPE algorithm
    fidelitytarg(H, psi0, psitarg, QI, QQ)


def fidelitytarg(H, psi0, psitarg, QI, QQ):
    # Schrodinger equation solver
    result = mesolve(H, psi0, times, args={'QI': QI, 'QQ': QQ})
    # result = mesolve()

    # Prints initial fidelity
    print("|<psi_0|psi_target>|² Initial fidelity: ",
          fidelity(psitarg, psitarg), "\n")

    # Bloch spheare with the initial, target and evaluated qubit states
    bl = Bloch()
    bl.add_states(ptrace(psi0, 0))  # Initial
    bl.add_states(ptrace(psitarg, 0))  # Target
    bl.add_states(ptrace(result.states[-1], 0))  # Evaluated
    bl.show()

    # Makes fidelity vs time plot
    fidel = np.zeros(len(result.states))
    for i in range(len(fidel)):
        fidel[i] = fidelity(result.states[i], psitarg)
    fig, ax = subplots()
    ax.plot(times, fidel)

    # Returns the final fidelity to the target
    return fidelity(result.states[-1], psitarg)


def setH0():
    H0 = 0.1*qd*q  # Detuning
    return H0


def qI(t, args):
    QI = args['QI']
    i = int(t/dt)
    return QI[min(np.abs(i), Ns-1)]


def qQ(t, args):
    QQ = args['QQ']
    i = int(t / dt)
    return QQ[min(np.abs(i), Ns-1)]


def settarg():
    return basis(2, 1)


if __name__ == '__main__':
    main()

print("\n---------------------------------------------------------------------\n")
