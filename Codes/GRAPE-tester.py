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


w = 5*2*np.pi*0
alpha = 0.2*np.pi*2
qubit_levels = 4
cavity_levels = 8
# Basic operators
q = np.array(destroy(qubit_levels))
qd = q.conj().T

c = np.array(destroy(cavity_levels))
cd = c.conj().T

q2 = q@q
qd2 = qd@qd

H0 = qd@q - (alpha/2)*qd2@q2
Hq_I = q+qd
Hq_Q = 1j*(q-qd)



# # Basic operators
# # Atom(a and a-dagger)
# q = tensor(destroy(2), qt.qeye(cavity_levels))
# qd = q.dag()
# # Cavity(c, c-dagger, and sigmaZ)
# c = tensor(qeye(2), qt.destroy(cavity_levels))
# cd = c.dag()

# # Hamiltonian operators
# H0 = qd * q  + cd*c + qd*q*cd*c# + (alpha/2)*qd2 @ q2
# Hq_I = q + qd
# Hq_Q = 1j * (q - qd)

# Hc_I = c + cd
# Hc_Q = 1j * (c - cd)

# H0 = np.array(H0)
# Hq_I = np.array(Hq_I)
# Hq_Q = np.array(Hq_Q)
# Hc_I = np.array(Hc_I)
# Hc_Q = np.array(Hc_Q)

# q = np.array(q)
# qd = np.array(qd)

# c = np.array(c)
# cd = np.array(cd)

print(H0.shape)
# Time variables
T = 10  # Total time of simulation
Ns = 5000 # Number of time steps
dt = T/Ns
times = np.linspace(0.0, T, Ns)

# initial and target qubit states
# qubit_initial = np.array(qt.basis(qubit_levels, 0))
# qubit_target = np.array(qt.basis(qubit_levels, 0))

# cavity_initial = np.array(qt.basis(cavity_levels, 0))
# cavity_target = np.array(qt.basis(cavity_levels, 6))

# psi_initial = np.kron(qubit_initial, cavity_initial)
# psi_target = np.kron(qubit_target, cavity_target)
# print(len(psi_initial))
#
# psi_initial = np.array([1,0])
# psi_target = np.array([0,1])
#
psi_initial = np.zeros(qubit_levels)
psi_target = np.zeros(qubit_levels)
psi_initial[0] = 1 
psi_target[1] = 1
# psi_initial = np.kron()
# psi_initial = np.reshape(psi_initial, (qubit_levels,1))
# psi_target = np.reshape(psi_target, (qubit_levels,1))

epsilon_max = 50
sigma = 3
# Create initial guess for the pulses
QI = np.sin(times)*(np.pi/(2*T)) + (np.random.random(Ns)-0.5)*0.1*np.exp(-((times -T/2)/sigma)**2)*0
QQ = np.cos(times)*(np.pi/(2*T)) + (np.random.random(Ns)-0.5)*0.1*np.exp(-((times -T/2)/sigma)**2)*0
# QI = (QI/np.max(QI))*epsilon_max*0.7*0 + epsilon_max/3
# QQ = (QQ/np.max(QQ))*epsilon_max*0.7*0 + epsilon_max/3

QI = np.sin(times)*(np.pi/(2*T))
QQ = np.cos(times)*(np.pi/(2*T))
# Creating a gaussian window for initial control amps guess
gaussian_window = gaussian(int(Ns/10), Ns/50, 1)
# Initial guess is a convolution of random number and gaussian window between
# The amplitude of the random vector in the convolution(from -rand_amp to rand_amp)
rand_amp_Q = 1/1000
rand_amp_I = 1/1000
# Calculating the convolutions themselves
conv_I = (ndi.convolve((np.random.random(Ns) - 0.5) * 2 * rand_amp_I, gaussian_window, mode='wrap'))
conv_Q = (ndi.convolve((np.random.random(Ns) - 0.5) * 2 * rand_amp_Q, gaussian_window, mode='wrap'))

QI = conv_I
QQ = conv_Q

guess_freq = 1
guess_width = 10
guess_amp = 1e-2
QI = guess_amp*np.exp((-times**2)/guess_width**2)*np.sin(times)
QQ = guess_amp*np.exp((-times**2)/guess_width**2)*np.cos(times)

pulse = np.array([QI, QQ])
pulse = pulse.flatten()

itime = time.time()
# Create the GrapePulse object :)
test_pulse = grape.GrapePulse(psi_initial, psi_target, T, Ns, H0, [Hq_I, Hq_Q], pulse, print_fidelity=True,
                              epsilon_max=epsilon_max, lambda_band_lin=0.002, lambda_amp_lin=0.0)
# test_pulse.cost_gradient(pulse, debug_fidelity=True)
# optimize with grape
pulse, fidelity = test_pulse.optimize()
print("Total time: ", time.time() - itime)

# Some graphs
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

# Plotting the final control pulses
axes[0, 1].set_title("final pulse")
axes[0, 1].step(times,  pulse[0])
axes[0, 1].step(times,  pulse[1])

FQI = np.fft.ifft(pulse[0])
FQQ = np.fft.ifft(pulse[1])
# Plotting the final control pulses in frequency space
axes[1, 1].set_title("final pulse frequency space")
axes[1, 1].step(fft_freq, FQI)
axes[1, 1].step(fft_freq, FQQ)

plt.show()
