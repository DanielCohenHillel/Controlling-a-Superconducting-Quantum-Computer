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
w = 1
alpha = 0.05/1000
qubit_levels = 3
epsilon_max = 50
sigma = 30

# -- Basic operators --
q = destroy(qubit_levels)
qd = q.dag()

q2 = q * q
qd2 = qd * qd

# -- Hamiltonians --
H0 = w*qd*q - (alpha/2)*qd2*q2
Hq_I = q + qd
Hq_Q = 1j*(q - qd)

# -- Time variables --
T = 250
Ns = 1000

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

# QI = conv_I
# QQ = conv_Q

pulse = np.array([QI, QQ])

itime = time.time()

# -- Using GrapePulse class --
# Create the GrapePulse object :)
test_pulse = grape.GrapePulse(psi_initial, psi_target, T, Ns, H0, [Hq_I, Hq_Q], pulse, constraints=True, print_fidelity=True,
                              max_amp=epsilon_max, lambda_band_lin=0.0, lambda_amp_lin=0, fix_amp_max=False)

# -- Test Gradient --
test_pulse.cost(pulse*2)
test_pulse.cost_gradient(pulse*2, debug_fidelity=True)

# -- Cost Function Optimization --
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
