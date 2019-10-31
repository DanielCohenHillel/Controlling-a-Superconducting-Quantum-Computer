# Imports
from qutip import *
from qutip.control import *
import qutip as qt
import numpy as np
from matplotlib import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots
import qutip.control.pulseoptim as cpo
from qutip.qip import hadamard_transform
import qutip.logging_utils as logging
import  scipy.optimize as spopt
from scipy.sparse.linalg import expm_multiply, expm

logger = logging.get_logger()
log_level = logging.INFO

print("\n---------------------------------------------------------------------\n")

h_bar = 1

# Basic operators
q = destroy(2)
qd = q.dag()

psi_initial = basis(2, 0)
psi_target = basis(2, 1)

# Time variables
Ns = 50  # Number of time steps
T = 5  # Total time of simulation
dt = T/Ns  # Time step size
times = np.linspace(0.0, T, Ns)  # All times array

# Constraints
epsilon_max = 0.3
epsilon_soft_max = 1


def main():
    H0 = 1.5*qd*q
    Hq_I = q + qd
    Hq_Q = 1j*(q-qd)

    H_d = H0
    H_c = [Hq_I, Hq_Q]

    HH = [H0, Hq_I, Hq_Q]

    U_0 = psi_initial
    U_targ = psi_target
    n_ts = Ns
    evo_time = T
    # Fidelity error target
    fid_err_targ = 1e-10
    # Maximum iterations for the optisation algorithm
    max_iter = 20000
    # Maximum (elapsed) time allowed in seconds
    max_wall_time = 120
    # Minimum gradient (sum of gradients squared)
    # as this tends to 0 -> local minima has been found
    min_grad = 1e-20
    p_type = 'ZERO'
    f_ext = "{}_n_ts{}_ptype{}.txt".format("Hi", n_ts, p_type)

    NN = 50

    result = qt.control.grape_unitary(U_targ, H_d, H_c, NN, times)

    QI = result.u[-1, 0, :] + np.random.random(len(result.u[-1, 0, :]))*0.02*0
    QQ = result.u[-1, 1, :] + np.random.random(len(result.u[-1, 0, :]))*0.02*0

    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.set_title("Qutip's grape")
    ax1.step(times, QI)
    ax1.step(times, QQ)

    # Initial guess is random
    QQ = (np.random.random(len(QQ))-1)*0 + 0.01
    QI = (np.random.random(len(QI))-1)*0 + 0.01

    control_vars = list(QI) + list(QQ)
    result = optimize_pulse_unitary_daniel(control_vars, HH, psi_initial, psi_target)

    a, b, c = result
    xI, xQ = np.split(np.array(a), 2)

    QI = epsilon_max * np.tanh(xI)
    QQ = epsilon_max * np.tanh(xQ)

    ax2.set_title("My Grape")
    ax2.step(times, QI)
    ax2.step(times, QQ)

    print("Final fidelity: ", (1-fidelitytarg(list(QI) + list(QQ), HH, psi_initial, psi_target))*100, "%")


def optimize_pulse_unitary_daniel(x0, H, psi_initial, psi_target):
    result = spopt.fmin_l_bfgs_b(fidelitytarg_constraints, x0, fidelitytarg_grad_constraints, maxiter=1000, pgtol=1e-10,
                                 args=(H, psi_initial, psi_target))
    return result


def fidelitytarg(control_pulses,  HH, psi_initial, psi_target):

    QI, QQ = np.split(np.array(control_pulses), 2)

    prod = 1
    for k in range(Ns):
        U_k = (-(1j * dt) * H(HH[0], HH[1], HH[2], QI, QQ, k)).expm()
        prod = U_k * prod
    psi_final = prod * psi_initial

    fid = qutip.fidelity(psi_target, psi_final)

    print("\n-  |<psi | psi_target>|²  fidelity: ", fid*100, "%", "\n")

    return 1-fid


def fidelitytarg_constraints(control_pulses,  HH, psi_initial, psi_target):

    xI, xQ = np.split(np.array(control_pulses), 2)

    QI = epsilon_max * np.tanh(xI)
    QQ = epsilon_max * np.tanh(xQ)

    prod = 1
    for k in range(Ns):
        U_k = (-(1j * dt) * H(HH[0], HH[1], HH[2], QI, QQ, k)).expm()
        prod = U_k * prod
    psi_final = prod * psi_initial

    fid = qutip.fidelity(psi_target, psi_final)

    print("\n-  |<psi | psi_target>|²  fidelity: ", fid*100, "%", "\n")

    return (1 - fid) - constraint(QI, QQ)


def fidelitytarg_grad(control_pulses, HH, psi_initial, psi_target):
    QI, QQ = np.split(np.array(control_pulses), 2)

    psi_bwd = []
    psi_fwd = []
    U_k = []

    for k in range(Ns):
        U_k.append((-1*(1j*dt/h_bar)*H(HH[0], HH[1], HH[2], QI, QQ, k)).expm())

    for k in range(Ns):
        if k == 0:
            psi_fwd.append(psi_initial)
            psi_bwd.append(psi_target)
        else:
            psi_fwd.append(U_k[k]*psi_fwd[k-1])
            psi_bwd.insert(0, (U_k[Ns-k].dag()*psi_bwd[0]))

    for k in range(-1, Ns):
        if k == -1:
            psi_bwd.append(psi_target)
        else:
            psi_bwd.insert(0, (U_k[Ns - k-1].dag() * psi_bwd[0]))

    c_final = []
    for k in range(2*Ns):
        c_final.append(np.abs((1j*dt/h_bar)*(psi_bwd[k % Ns].dag()*HH[1 if k < Ns else 2]*psi_fwd[k % Ns]).data[0, 0]))
    return np.transpose(c_final)


def fidelitytarg_grad_constraints(control_pulses, HH, psi_initial, psi_target):
    xI, xQ = np.split(np.array(control_pulses), 2)

    QI = epsilon_max * np.tanh(xI)
    QQ = epsilon_max * np.tanh(xQ)

    psi_bwd = []
    psi_fwd = []
    U_k = []

    for k in range(Ns):
        U_k.append((-1*(1j*dt/h_bar)*H(HH[0], HH[1], HH[2], QI, QQ, k)).expm())

    prod = 1
    for k in range(Ns):
        Uk = ((-1j * dt / h_bar) * H(HH[0], HH[1], HH[2], QI, QQ, k)).expm()
        prod = Uk * prod
    psi_final = prod * psi_initial
    # print(prod*psi_initial)
    # print(fidelitytarg(list(QI) + list(QQ), HH, psi_initial, psi_target))

    for k in range(Ns):
        if k == 0:
            psi_fwd.append(psi_initial)
            psi_bwd.append(psi_target)
        else:
            psi_fwd.append(U_k[k]*psi_fwd[k-1])
            psi_bwd.insert(0, (U_k[Ns-k].dag()*psi_bwd[0]))

    for k in range(-1, Ns):
        if k == -1:
            psi_bwd.append(psi_target)
        else:
            psi_bwd.insert(0, (U_k[Ns - k -1].dag() * psi_bwd[0]))

    c_final = []
    for k in range(2*Ns):
        c_final.append(np.abs((-1j*dt/h_bar)*(psi_bwd[k % Ns].dag()*HH[1 if k < Ns else 2]*psi_fwd[k % Ns]).data[0, 0])
                       / (np.cosh((xI[k] if k < Ns else xQ[k-Ns])**2)))
    return epsilon_max*fidelity(psi_target, psi_final)*np.transpose(c_final)


def constraint(QI, QQ):
    constraint_total = 0
    lambda_amp = 0.0001
    lambda_band = -1e10
    for i in range(len(QI)):
        # Amp constraint
        # constraint_total = constraint_total + (QI[i]**2 + QQ[i]**2)*lambda_amp

        # Bandwith constraint
        # if i != (len(QI)-1):
        #     constraint_total = constraint_total + (np.exp(((QI[i + 1] - QI[i]) ** 2 + (QQ[i + 1] - QQ[i]) ** 2)
        #                                                   / (epsilon_soft_max**2))-1)*lambda_band
        if i != (len(QI) - 1):
            constraint_total = constraint_total + lambda_band*(np.abs((QI[i + 1] - QI[i]) ** 2
                                                                     + (QQ[i + 1] - QQ[i]) ** 2))
    print("Constraint Total: ", constraint_total)
    return constraint_total


def H(H0, Hq_I, Hq_Q, QI, QQ, k):
    return H0 + Hq_I*QI[k] + Hq_Q*QQ[k]


def qI(t, args):
    QI = args['QI']
    i = int(t/dt)
    return QI[min(np.abs(i), Ns-1)]


def qQ(t, args):
    QQ = args['QQ']
    i = int(t / dt)
    return QQ[min(np.abs(i), Ns-1)]


if __name__ == '__main__':
    main()

print("\n---------------------------------------------------------------------\n")

# fig1 = plt.figure()
# ax1 = fig1.add_subplot(2, 1, 1)
# ax1.set_title("Initial control amps")
# # ax1.set_xlabel("Time")
# ax1.set_ylabel("Control amplitude")
# ax1.step(result.time,
#          np.hstack((result.initial_amps[:, 0], result.initial_amps[-1, 0])),
#          where='post')
# ax1.step(result.time,
#          np.hstack((result.initial_amps[:, 1], result.initial_amps[-1, 0])),
#          where='post')
#
# ax2 = fig1.add_subplot(2, 1, 2)
# ax2.set_title("Optimised Control Sequences")
# ax2.set_xlabel("Time")
# ax2.set_ylabel("Control amplitude")
# ax2.step(result.time,
#          np.hstack((result.final_amps[:, 0], result.final_amps[-1, 0])),
#          where='post')
# ax2.step(result.time,
#          np.hstack((result.final_amps[:, 1], result.final_amps[-1, 0])),
#          where='post')
#
# plt.tight_layout()
# plt.show()
