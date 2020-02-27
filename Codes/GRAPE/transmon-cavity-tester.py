from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import grape  # This is mine :)
import time
import scipy.ndimage as ndi
import scipy.signal

plt.style.use("fivethirtyeight")


def smooth(y, box_pts):
    box = np.ones(int(box_pts))/int(box_pts)
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


# Time variables
T = 10  # Total time of simulation
Ns = 100  # Number of time steps
dt = T/Ns
times = np.linspace(0.0, T, Ns)

# Values from Phillip's thesis
w_a = 5.66*0.5
w_c = 4.5
chi = 2

qubit_levels = 2
cavity_levels = 10

psi_initial = tensor(basis(qubit_levels, 0), basis(cavity_levels, 0))
psi_target = tensor(basis(qubit_levels, 0), basis(cavity_levels, 1))

# Anharmonic
alpha = w_a*0.05
K = w_c*10**(-3)
chitag = chi*10**(-2)

guess_freq = 1
guess_width = 50
guess_amp = (np.pi/(2*T))*2
QI = guess_amp*np.exp((-(times)**2)/guess_width**2) * \
    np.sin(w_a*times)*(1 + 0.2*np.random.random(Ns))
QQ = guess_amp*np.exp((-(times)**2)/guess_width**2) * \
    np.cos(w_a*times)*(1 + 0.2*np.random.random(Ns))

CI = guess_amp*np.exp((-(times)**2)/guess_width**2) * \
    np.sin(w_c*times)*(1 + 0.2*np.random.random(Ns))
CQ = guess_amp*np.exp((-(times)**2)/guess_width**2) * \
    np.cos(w_c*times)*(1 + 0.2*np.random.random(Ns))


# Basic operators
a = tensor(destroy(qubit_levels), qeye(cavity_levels))
ad = a.dag()

c = tensor(qeye(qubit_levels), destroy(cavity_levels))
cd = c.dag()

Ha = w_a*ad*a + (alpha/2)*ad*ad*a*a
Hc = w_c*cd*c + (K/2)*cd*cd*c*c
Hi = chi*cd*c*ad*a + chitag*cd*cd*c*c*ad*a
H0 = Ha + Hc + Hi

Ha_I = a + ad
Ha_Q = 1j*(a - ad)

Hc_I = c + cd
Hc_Q = 1j*(c - cd)

drive_hamiltonians = [Ha_I, Ha_Q, Hc_I, Hc_Q]

epsilon_max = 50

pulse = np.array([QI, QQ, CI, CQ])
pulse = pulse.flatten()
itime = time.time()
# Create the GrapePulse object :)
test_pulse = grape.GrapePulse(psi_initial, psi_target, T, Ns, H0, drive_hamiltonians, pulse, print_fidelity=True,
                              max_amp=epsilon_max, lambda_band_lin=0.0, lambda_amp_lin=0.0, fix_amp_max=False)
# optimize with grape
pulse, fidelity = test_pulse.optimize()
print(np.abs(fidelity))
# while np.abs(fidelity) < 0.99:
#     print("Bad Fidelity, Retrying")
#     r = np.random.random(Ns*len(drive_hamiltonians))
#     # r = scipy.signal.savgol_filter(r - 0.5, int((Ns/20) + 1-(Ns/20)%2), 3)
#     r = r-0.5  # smooth(r-0.5, Ns/10)
#     test_pulse.initial_pulse = pulse.flatten() + r*(2*np.max(pulse)*0.8)
#     pulse, fidelity = test_pulse.optimize()
print("Total time: ", time.time() - itime)

# # Some graphs
# # Creating plots for the amplitudes
# fig, axes = plt.subplots(2, 2)
# # print(axes)
# axes[0, 0].set_title("Initial pulse")

# axes[0, 0].step(times, QI)
# axes[0, 0].step(times, QQ)

# FQI = np.fft.ifft(QI)
# FQQ = np.fft.ifft(QQ)

# fft_freq = np.fft.fftfreq(Ns, dt)
# # Plotting the final control pulses in frequency space
# axes[1, 0].set_title("final pulse frequency space")
# axes[1, 0].step(fft_freq, FQI)
# axes[1, 0].step(fft_freq, FQQ)

# # Plotting the final control pulses
# axes[0, 1].set_title("final pulse")
# axes[0, 1].step(times,  pulse[0])
# axes[0, 1].step(times,  pulse[1])

# FQI = np.fft.ifft(pulse[0])
# FQQ = np.fft.ifft(pulse[1])
# # Plotting the final control pulses in frequency space
# axes[1, 1].set_title("final pulse frequency space")
# axes[1, 1].step(fft_freq, FQI)
# axes[1, 1].step(fft_freq, FQQ)

# # plt.show()
# # --- Cavity ---
# fig, axes_cav = plt.subplots(2)
# axes_cav[0].set_title("Initial Cavity Pulses")
# # axes_cav[0].legend(["I", "Q"])
# axes_cav[0].plot(times, CI)
# axes_cav[0].plot(times, CQ)

# axes_cav[1].set_title("Final Cavity Pulses")
# axes_cav[0].legend([r'$\epsilon_I^{Cavity}$', r'$\epsilon_Q^{Cavity}$'])
# axes_cav[1].plot(times, pulse[2])
# axes_cav[1].plot(times, pulse[3])

# # --- Transmon ---
# fig_cav, axes = plt.subplots(2)
# axes[0].set_title("Initial Qubit Pulses")
# axes[0].plot(times, QI)
# axes[0].plot(times, QQ)
# axes[0].legend([r'$\epsilon_I^{Qubit}$', r'$\epsilon_Q^{Qubit}$'])
# # axes[0].legend(['QI', 'QQ'])

# # -- Final --
# axes[1].set_title("Final Qubit Pulses")
# axes[1].plot(times,  pulse[0])
# axes[1].plot(times,  pulse[1])


# axes[0].set_xlabel("Time (ns)")
# axes[1].set_xlabel("Time (ns)")
# axes_cav[0].set_xlabel("Time (ns)")
# axes_cav[1].set_xlabel("Time (ns)")

# axes[0].set_ylabel("Amplitude (MHz)")
# axes[1].set_ylabel("Amplitude (MHz)")
# axes_cav[0].set_ylabel("Amplitude (MHz)")
# axes_cav[1].set_ylabel("Amplitude (MHz)")

# # plt.tight_layout()
# fig.tight_layout()
# fig_cav.tight_layout()
plt.show()
