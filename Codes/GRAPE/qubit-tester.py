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


w = 1
alpha = 0.05
qubit_levels = 3

# Basic operators
q = destroy(qubit_levels)
qd = q.dag()

q2 = q * q
qd2 = qd * qd

H0 = w*qd*q - (alpha/2)*qd2*q2
Hq_I = q + qd
Hq_Q = 1j*(q - qd)

# Time variables
T = 1300  # Total time of simulation
Ns = 200  # Number of time steps
dt = T/Ns
times = np.linspace(0.0, T, Ns)

psi_initial = qt.basis(qubit_levels, 0)
psi_target = qt.basis(qubit_levels, 1)

epsilon_max = 50
sigma = 30

# Create initial guess for the pulses
QI = np.sin(times)*(np.pi/(2*T)) + (np.random.random(Ns)-0.5) * \
    0.1*np.exp(-((times - T/2)/sigma)**2)*0
QQ = np.cos(times)*(np.pi/(2*T)) + (np.random.random(Ns)-0.5) * \
    0.1*np.exp(-((times - T/2)/sigma)**2)*0

QI = np.sin(w*times)*(np.pi/(2*T))*1
QQ = np.cos(w*times)*(np.pi/(2*T))*1

# Creating a gaussian window for initial control amps guess
gaussian_window = gaussian(int(Ns/10), Ns/50, 1)

# Initial guess is a convolution of random number and gaussian window between
# The amplitude of the random vector in the convolution(from -rand_amp to rand_amp)
rand_amp_Q = 1/10
rand_amp_I = 1/10
# Calculating the convolutions themselves
conv_I = (ndi.convolve((np.random.random(Ns) - 0.5) *
                       2 * rand_amp_I, gaussian_window, mode='wrap'))
conv_Q = (ndi.convolve((np.random.random(Ns) - 0.5) *
                       2 * rand_amp_Q, gaussian_window, mode='wrap'))

QI = conv_I
QQ = conv_Q

guess_freq = 1
guess_width = T/2
guess_amp = (np.pi/(2*T))*2
# QI = guess_amp*np.exp((-times**2)/guess_width**2)*np.sin(w*times) + (np.random.random(Ns) - 0.5)*0.4
# QQ = guess_amp*np.exp((-times**2)/guess_width**2)*np.cos(w*times) + (np.random.random(Ns) - 0.5)*0.4

pulse = np.array([QI, QQ])
# pulse = pulse.flatten()

itime = time.time()
# Create the GrapePulse object :)
test_pulse = grape.GrapePulse(psi_initial, psi_target, T, Ns, H0, [Hq_I, Hq_Q], pulse, constraints=True, print_fidelity=True,
                              max_amp=epsilon_max, lambda_band_lin=3.1, lambda_amp_lin=1*0, fix_amp_max=False)
# test_pulse.cost(pulse*2)
# test_pulse.cost_gradient(pulse*2, debug_fidelity=True)
# optimize with grape
pulse, fidelity = test_pulse.optimize()
print("Total time: ", time.time() - itime)
# test_pulse.run_operator(pulse)
# Some graphs
# Creating plots for the amplitudes
fig, axes = plt.subplots(2)
# print(axes)
axes[0].set_title("Initial pulse")
axes[0].plot(times, QI)
axes[0].plot(times, QQ)

# FQI = np.fft.ifft(QI)
# FQQ = np.fft.ifft(QQ)

# fft_freq = np.fft.fftfreq(Ns, dt)
# # Plotting the final control pulses in frequency space
# axes[1, 0].set_title("final pulse frequency space")
# axes[1, 0].step(fft_freq, FQI)
# axes[1, 0].step(fft_freq, FQQ)

# Plotting the final control pulses
axes[1].set_title("final pulse")
axes[1].plot(times,  pulse[0])
axes[1].plot(times,  pulse[1])

# FQI = np.fft.ifft(pulse[0])
# FQQ = np.fft.ifft(pulse[1])
# # Plotting the final control pulses in frequency space
# axes[1, 1].set_title("final pulse frequency space")
# axes[1, 1].step(fft_freq, FQI)
# axes[1, 1].step(fft_freq, FQQ)

plt.show()
