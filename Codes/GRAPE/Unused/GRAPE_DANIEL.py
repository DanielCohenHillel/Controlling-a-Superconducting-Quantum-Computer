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
import scipy.optimize as spopt
from scipy.sparse.linalg import expm_multiply, expm
import scipy.ndimage as ndi

print("\n---------------------------------------------------------------------\n")

# Basic operators
q = destroy(2)
qd = q.dag()

psi_initial = (basis(2, 0) + basis(2, 1))/np.sqrt(2)
# psi_initial = basis(2, 0)
psi_target = basis(2, 1)

# Time variables
Ns = 500  # Number of time steps
T = 10  # Total time of simulation
dt = T / Ns  # Time step size
times = np.linspace(0.0, T, Ns)  # All times array

# Constraints
epsilon_max = 0.5
epsilon_soft_max = 1

lambda_band_lin = 0.1

# Hamiltonian operators
H0 = qd * q  # Base Hamiltonian
Hq_I = q + qd
Hq_Q = 1j * (q - qd)
H_c = [Hq_I, Hq_Q]  # Drive Hamiltonians in list form for qutip GRAPE
HH = [H0, Hq_I, Hq_Q]  # All Hamiltonian components in list form

A = 0.5


def main():
    # Creating a gaussian window for initial control amps guess
    gaussian_window = gaussian(min(80, Ns), 40, 1, True)
    # Initial guess is a convolution of random number and gaussian window between
    # The amplitude of the random vector in the convolution(from -rand_amp to rand_amp)
    rand_amp_Q = .5
    rand_amp_I = .5
    # Calculating the convolutions themselves
    conv_I = (ndi.convolve((np.random.random(Ns) - 0.5) *
                           2 * rand_amp_I, gaussian_window, mode='wrap'))
    conv_Q = (ndi.convolve((np.random.random(Ns) - 0.5) *
                           2 * rand_amp_Q, gaussian_window, mode='wrap'))

    gauss_width = 0.1
    gauss = np.exp(-gauss_width*((times-T/2)**2)) - \
        np.exp(-gauss_width * ((T/2)**2))
    plt.figure()
    plt.plot(times, gauss)
    # Adding the convolutions to sin/cos waves
    QI = (np.sin(times-np.pi/4)*(np.pi/(2*T)) +
          (np.random.random(Ns)-0.5)*0.03)*gauss
    QQ = (np.cos(times+np.pi/4)*(np.pi/(2*T)) +
          (np.random.random(Ns)-0.5)*0.03)*gauss
    print(run_operator(QI, QQ))

    # QI = (np.random.random(Ns)-0.5)
    # QQ = (np.random.random(Ns)-0.5)

    # Creating plots for the amplitudes
    fig, axes = plt.subplots(2, 2)
    # print(axes)
    axes[0, 0].set_title("Initial pulse")
    axes[0, 0].step(times, QI)
    axes[0, 0].step(times, QQ)

    FQI = np.fft.ifft(QI)
    FQQ = np.fft.ifft(QQ)

    fft_freq = np.fft.fftfreq(Ns, dt)
    # Plotting the final control pulses in frequency space
    axes[1, 0].set_title("final pulse frequency space")
    axes[1, 0].step(fft_freq, FQI)
    axes[1, 0].step(fft_freq, FQQ)

    # Putting the two control pulses into one long array to send to the optimization algorithm
    control_vars = np.concatenate((QI, QQ))
    # Runs the optimization algorithm
    result, final_fidelity = optimize_pulse_unitary_daniel(control_vars)
    print("Final Fidelity = ", final_fidelity)
    # for i in range(5):  # TODO: Change 5 to a variable
    #     if np.abs(final_fidelity) > 0.99:
    #         print("Nice fidelity :)")
    #         break
    #     print(" ----- Trying again :(")
    #     result = optimize_pulse_unitary_daniel(control_vars + (np.random.random(len(control_vars)) - 0.5)*0.3)

    # Getting each pulse separately from the results
    xI, xQ = np.split(result, 2)

    # The actual control amplitudes after the amp hard-cutoff
    QI = epsilon_max * np.tanh(xI)
    QQ = epsilon_max * np.tanh(xQ)

    # Plotting the final control pulses
    axes[0, 1].set_title("final pulse")
    axes[0, 1].step(times, QI)
    axes[0, 1].step(times, QQ)

    FQI = np.fft.ifft(xI)
    FQQ = np.fft.ifft(xQ)
    # Plotting the final control pulses in frequency space
    axes[1, 1].set_title("final pulse frequency space")
    axes[1, 1].step(fft_freq, FQI)
    axes[1, 1].step(fft_freq, FQQ)

    run_operator(QI, QQ, show_bloch=True)
    print("\n\n---> Final fidelity: ",
          str(-1*get_fidelity(np.concatenate((QI, QQ)))*100) + "%")


def optimize_pulse_unitary_daniel(x0, constraints=True, fix_amp_max=True):
    # Using the L-BFGS-B optimization algorithm to find the minimum of the cost function
    if constraints:
        # With Constraints
        if fix_amp_max:
            result = spopt.fmin_l_bfgs_b(get_fidelity_constraints,
                                         np.arctanh(x0/epsilon_max), get_fidelity_gradient_constraints, factr=1e12)
        else:
            result = spopt.fmin_l_bfgs_b(
                get_fidelity_constraints, x0, get_fidelity_gradient_constraints)
    else:
        # Without Constraints
        result = spopt.fmin_l_bfgs_b(get_fidelity, x0, get_fidelity_gradient)
    return result[0:2]


# ----------------------- Without Constraints -----------------------------

def get_fidelity(control_pulses):
    # Get each control pulse separately
    QI, QQ = np.split(control_pulses, 2)

    # Calculates the final state and the fidelity from the run_operator function
    psi_final, fid = run_operator(QI, QQ)

    return -fid


def get_fidelity_gradient(control_pulses, debug_fidelity=False):
    # Get each control pulse separately
    QI, QQ = np.split(control_pulses, 2)

    # Initialize the elements used to calculate the gradient efficiently
    psi_bwd = []
    psi_fwd = []
    U_k = []

    # print("--------------------> ", H(H0, QI, QQ, 1).expm())
    # print("--------------------> ", QI[2])
    for k in range(Ns):
        U_k.append((1j * dt * H(H0, QI, QQ, k)).expm())

    for k in range(Ns):
        if k == 0:
            psi_fwd.append(psi_initial)
            # psi_bwd.append(psi_target)
        else:
            psi_fwd.append(U_k[k-1] * psi_fwd[k - 1])
            # psi_bwd.insert(0, (U_k[Ns - k].dag() * psi_bwd[0]))

    psi_bwd.append(1)
    for k in range(Ns):
        psi_bwd.insert(0, psi_bwd[0] * U_k[Ns - 1 - k])
    for k in range(Ns + 1):
        psi_bwd[k] = psi_target.dag() * psi_bwd[k]

    prod = 1
    for k in range(Ns):
        U_k = ((1j * dt) * H(H0, QI, QQ, k)).expm()
        prod = U_k * prod
    c = psi_target.dag() * prod * psi_initial

    c_final = []
    for k in range(Ns):
        c_final.append(2 * np.real(c * np.conjugate(1j * dt *
                                                    (psi_bwd[k] * Hq_I * psi_fwd[k])[0][0])))
    for k in range(Ns):
        c_final.append(2 * np.real(c * np.conjugate(1j * dt *
                                                    (psi_bwd[k] * Hq_Q * psi_fwd[k])[0][0])))

    if debug_fidelity:
        fig, axes = subplots(2, 2)
        axes[0, 0].set_title("gradient")
        axes[0, 0].plot(times, c_final[0:Ns])
        axes[0, 0].plot(times, c_final[Ns:2*Ns])
        axes[0, 0].plot(times, np.zeros(Ns))

        axes[1, 0].set_title("rough estimation of gradient")
        grad = - \
            spopt.approx_fprime(np.concatenate((QI, QQ)), get_fidelity, 0.001)
        axes[1, 0].plot(times, grad[0:Ns])
        axes[1, 0].plot(times, grad[Ns:])
        axes[1, 0].plot(times, np.zeros(Ns))

        axes[0, 1].set_title("QI")
        axes[0, 1].plot(times, QI)
        axes[0, 1].plot(times, QI/A, '--')
        axes[0, 1].plot(times, np.zeros(Ns))

        axes[1, 1].set_title("QQ")
        axes[1, 1].plot(times, QQ)
        axes[1, 1].plot(times, QQ / A, '--')
        axes[1, 1].plot(times, np.zeros(Ns))

    return np.transpose(c_final)


# ----------------------- With Constraints -----------------------------

def get_fidelity_constraints(control_pulses):
    # Get each control pulse separately
    xI, xQ = np.split(np.array(control_pulses), 2)
    xI[0] = 0
    xI[Ns-1] = 0
    xQ[0] = 0
    xQ[Ns-1] = 0

    # Hard cut-off amplitude for the control pulses
    QI = epsilon_max * np.tanh(xI)
    QQ = epsilon_max * np.tanh(xQ)

    # Calculates the final state and the fidelity from the run_operator function
    psi_final, fid = run_operator(QI, QQ)

    # Adding the 'soft' constraints
    fid_constraints = -fid + constraint(xI, xQ)
    return fid_constraints


def get_fidelity_gradient_constraints(control_pulses, debug_fidelity=False):
    # Get each control pulse separately
    xI, xQ = np.split(control_pulses, 2)
    xI[0] = 0
    xI[Ns - 1] = 0
    xQ[0] = 0
    xQ[Ns - 1] = 0
    # Hard cut-off amplitude for the control pulses
    QI = epsilon_max * np.tanh(xI)
    QQ = epsilon_max * np.tanh(xQ)

    # Initialize the elements used to calculate the gradient efficiently
    psi_bwd = []
    psi_fwd = []
    U_k = []

    # print("--------------------> ", H(H0, QI, QQ, 1).expm())
    # print("--------------------> ", QI[2])
    for k in range(Ns):
        U_k.append((1j * dt * H(H0, QI, QQ, k)).expm())

    for k in range(Ns):
        if k == 0:
            psi_fwd.append(psi_initial)
            # psi_bwd.append(psi_target)
        else:
            psi_fwd.append(U_k[k - 1] * psi_fwd[k - 1])
            # psi_bwd.insert(0, (U_k[Ns - k].dag() * psi_bwd[0]))

    psi_bwd.append(1)
    for k in range(Ns):
        psi_bwd.insert(0, psi_bwd[0] * U_k[Ns - 1 - k])
    for k in range(Ns + 1):
        psi_bwd[k] = psi_target.dag() * psi_bwd[k]

    prod = 1
    for k in range(Ns):
        U_k = ((1j * dt) * H(H0, QI, QQ, k)).expm()
        prod = U_k * prod
    c = psi_target.dag() * prod * psi_initial

    c_final = []
    for k in range(Ns):
        c_final.append(2 * np.real(c * np.conjugate(1j * dt *
                                                    (psi_bwd[k] * Hq_I * psi_fwd[k])[0][0])))
    for k in range(Ns):
        c_final.append(2 * np.real(c * np.conjugate(1j * dt *
                                                    (psi_bwd[k] * Hq_Q * psi_fwd[k])[0][0])))

    xIxQ = np.concatenate((xI, xQ))
    gradient = (np.asarray(np.transpose(c_final)) -
                constraint_gradient(xI, xQ))*(epsilon_max / (np.cosh(xIxQ)**2))
    if debug_fidelity:
        fig, axes = subplots(2, 2)
        axes[0, 0].set_title("gradient")
        axes[0, 0].plot(times, np.transpose(gradient)[0:Ns])
        axes[0, 0].plot(times, np.transpose(gradient)[Ns:2 * Ns])
        axes[0, 0].plot(times, np.zeros(Ns))

        axes[1, 0].set_title("rough estimation of gradient")
        grad = -spopt.approx_fprime(np.concatenate((xI, xQ)),
                                    get_fidelity_constraints, 0.00001)
        axes[1, 0].plot(times, grad[0:Ns])
        axes[1, 0].plot(times, grad[Ns:])
        axes[1, 0].plot(times, np.zeros(Ns))

        axes[0, 1].set_title("QI")
        axes[0, 1].plot(times, QI)
        axes[0, 1].plot(times, QI / A, '--')
        axes[0, 1].plot(times, np.zeros(Ns))

        axes[1, 1].set_title("QQ")
        axes[1, 1].plot(times, QQ)
        axes[1, 1].plot(times, QQ / A, '--')
        axes[1, 1].plot(times, np.zeros(Ns))

    return -gradient


def constraint(QI, QQ):
    constraint_total = 0
    lambda_amp = 0
    for k in range(Ns-1):
        constraint_total = constraint_total + \
            lambda_band_lin*(QI[k+1] - QI[k])**2
        constraint_total = constraint_total + \
            lambda_band_lin*(QQ[k+1] - QQ[k])**2
    print("\n-)Constraint Total: ", constraint_total)
    return constraint_total


def constraint_gradient(QI, QQ):
    constraint_total = 0
    lambda_amp = 0
    g_band_lin = np.zeros(2*Ns)
    g_band_lin[0] = -2*(QI[1] - QI[0])
    g_band_lin[Ns - 1] = 2*(QI[Ns-1] - QI[Ns-2])
    g_band_lin[0] = -2 * (QQ[1] - QQ[0])
    g_band_lin[Ns - 1] = 2 * (QQ[Ns - 1] - QQ[Ns - 2])
    for k in range(1, Ns-1):
        g_band_lin[k] = 4*QI[k] - 2*(QI[k+1] + QI[k-1])
    for k in range(1, Ns-1):
        g_band_lin[k + Ns] = 4*QQ[k] - 2*(QQ[k+1] + QQ[k-1])
    # for i in range(Ns - 1):
    #     if i != (Ns - 1):
    #         constraint_total = constraint_total \
    #                            + lambda_band_lin * (np.abs((QI[i + 1] - QI[i]) ** 2 + (QQ[i + 1] - QQ[i]) ** 2)) \
    #                            + lambda_amp * (QI[i] ** 2 + QQ[i] ** 2)
    constraint_total = lambda_band_lin*g_band_lin
    # print("\n-)Constraint Gradient Total: ", constraint_total)
    return constraint_total


def run_operator(QI, QQ, show_bloch=False, calc_fidelity=True):
    prod = 1
    if show_bloch:
        b = Bloch()
        b.add_states(psi_initial)
        for k in range(Ns):
            U_k = ((1j * dt) * H(H0, QI, QQ, k)).expm()
            prod = U_k * prod
            b.add_states(prod * psi_initial, kind='point')
        psi_final = prod * psi_initial
        b.add_states(psi_final)
        b.show()
    else:
        for k in range(Ns):
            U_k = ((1j * dt) * H(H0, QI, QQ, k)).expm()
            prod = U_k * prod
        psi_final = prod * psi_initial

    if calc_fidelity:
        fid = (np.abs(psi_target.dag() * psi_final)**2)[0]
        print("\n-> Fidelity: ", fid)
        return psi_final, fid
    return psi_final


def run_operatorr(QI, QQ, calc_fidelity=True):
    C = []
    prod = 1
    for k in range(Ns):
        U_k = ((1j * dt) * H(H0, QI, QQ, k)).expm()
        prod = U_k * prod
        C.append(np.abs((prod*psi_initial)[1][0]))
    psi_final = prod*psi_initial

    print(np.argmax(C))
    plt.title("|1>")
    plt.plot(times, C[:])

    if calc_fidelity:
        fid = fidelity(psi_final, psi_target)
        print("\n-> Fidelity: ", fid)
        return psi_final, fid
    return psi_final


def H(H0, QI, QQ, k):
    return H0 + Hq_I * QI[k] + Hq_Q * QQ[k]


def gaussian(size, sigma, amp, graph=False):
    gaussian_window = np.zeros(size)
    for x in range(size):
        gaussian_window[x] = amp * np.exp(-(x - size / 2 ** 2) / sigma ** 2)
    if graph:
        plt.figure()
        plt.plot(times[0:size], gaussian_window)
    return gaussian_window


if __name__ == '__main__':
    main()
    # guess = np.concatenate((np.sin(times)*A + 0*(np.random.random(Ns)-0.5), np.cos(times)*A + 0*(np.random.random(Ns)-0.5)))
    # get_fidelity_gradient_constraints(guess, debug_fidelity=True)

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
#
# def qI(t, args):
#     QI = args['QI']
#     i = int(t / dt)
#     return QI[min(np.abs(i), Ns - 1)]
#
#
# def qQ(t, args):
#     QQ = args['QQ']
#     i = int(t / dt)
#     return QQ[min(np.abs(i), Ns - 1)]
