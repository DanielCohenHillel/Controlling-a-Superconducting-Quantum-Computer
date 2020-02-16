from qutip import *
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import grape  # This is mine :)
import time
import scipy.ndimage as ndi


def gaussian(size, sigma, amp, graph=False):
    gaussian_window = np.zeros(size)
    for x in range(size):
        gaussian_window[x] = amp * np.exp(-(x - size / 2 ** 2) / sigma ** 2)
        if graph:
            plt.figure()
            plt.plot(times[0:size], gaussian_window)
    return gaussian_window


# -- Constants --
w = 2*np.pi*0
alpha = 2*np.pi*0.15
epsilon_max = 50
sigma = 30
qubit_levels = 2

# -- Basic operators --
q = destroy(qubit_levels)
qd = q.dag()

q2 = q * q
qd2 = qd * qd

# -- Hamiltonians --
H0 = w*qd*q - (alpha/2)*qd2*q2
Hq_I = (q + qd)*0.5
Hq_Q = 0.5*1j*(q - qd)

# -- Time variables --
T = 30
Ns = 200

dt = T/Ns
times = np.linspace(0.0, T, Ns)

# -- Initial and Final States --
psi_initial = qt.basis(qubit_levels, 0)
psi_target = qt.basis(qubit_levels, 1)


# -- Initial Pulses --
# -- sin\cos Solution Pulses --
QI = np.sin(w*times)*(np.pi/(2*T))
QQ = np.cos(w*times)*(np.pi/(2*T))

# -- Random Initial Pulses --
gaussian_window = gaussian(int(Ns/10), Ns/50, 1)

rand_amp_Q = 1/10
rand_amp_I = 1/10

conv_I = (ndi.convolve((np.random.random(Ns) - 0.5) *
                       2 * rand_amp_I, gaussian_window, mode='wrap'))
conv_Q = (ndi.convolve((np.random.random(Ns) - 0.5) *
                       2 * rand_amp_Q, gaussian_window, mode='wrap'))

QI = conv_I
QQ = conv_Q*0.03

pulse = np.array([QI, QQ])

itime = time.time()

# -- Using GrapePulse class --
# Create the GrapePulse object :)
test_pulse = grape.GrapePulse(psi_initial, psi_target, T, Ns, H0, [Hq_I, Hq_Q], pulse, constraints=True, print_fidelity=True,
                              max_amp=epsilon_max, lambda_band_lin=0.0, lambda_amp_lin=0, fix_amp_max=False)

# -- Test Gradient --
# test_pulse.cost(pulse*2)
# test_pulse.cost_gradient(pulse*2, debug_fidelity=True)

# -- Cost Function Optimization --
# pulse, fidelity = test_pulse.optimize()

print("Going to DRAG")
qubit_levels = 3

# -- Basic operators --
q = destroy(qubit_levels)
qd = q.dag()

q2 = q * q
qd2 = qd * qd

# -- Hamiltonians --
H0 = w*qd*q - (alpha/2)*qd2*q2
Hq_I = 0.5*(q + qd)
Hq_Q = 0.5*1j*(q - qd)

# -- Initial and Final States --
psi_initial = qt.basis(qubit_levels, 0)
psi_target = qt.basis(qubit_levels, 1)

# QI = pulse[0]
# QQ = pulse[1]
QI = np.sin(w*times)
# QQ = np.cos(w*times)
sig = 4.0
A = np.sqrt(np.pi/2)/(sig)
QI = np.exp((-(times-T/2)**2)/(2*(sig)**2))*A*0 + np.sin(times)
QQ = (-0.5*np.exp((-(times-T/2)**2)/(2*(sig)**2)) *
      A*(times-T/2)/(alpha*sig**2))*0 + np.cos(times)

#   0*np.exp((-(T/2)**2)/(2*(sig)**2)))*A
pulse = np.array([QI, QQ])*2


test_pulse = grape.GrapePulse(psi_initial, psi_target, T, Ns, H0, [Hq_I, Hq_Q], pulse, constraints=False, print_fidelity=True,
                              max_amp=epsilon_max, lambda_band_lin=0.0, lambda_amp_lin=0, fix_amp_max=False, drag=True)

# -- Cost Function Optimization --
print("3-level before DRAG: ", test_pulse.run_operator(pulse, show_prob=True))
print("3-level before DRAG: ", test_pulse.run_operator(pulse))
test_pulse.cost(pulse*2)
test_pulse.cost_gradient(pulse*2, debug_fidelity=True)
# pulse, fidelity = test_pulse.optimize()

print("Total time: ", time.time() - itime)

# --- Pulses Graphs ---
# -- Initial --
fig, axes = plt.subplots(2)
axes[0].set_title("Initial pulse")
axes[0].plot(times, QI)
axes[0].plot(times, QQ)
# -- Final --
axes[1].set_title("final pulse")
axes[1].plot(times,  pulse[0])
axes[1].plot(times,  pulse[1])

plt.show()
