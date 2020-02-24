import numpy as np
import qutip as qt
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')

# N = 1000
# W = 50
# x = np.linspace(-1, 6, N)
# w = 3
# sig = 1
# gauss = np.exp((-(x-w)**2)/0.3)

T = 10
Nt = 1000
t = np.linspace(0, T, Nt)
# fourier = np.exp((-t**2)/4 - 1j*w*t)

# plt.plot(t, np.real(fourier))
# plt.plot(t, np.imag(fourier))
# plt.plot(t, np.exp((-t**2)/4), '--')
# R = np.real(fourier)
# I = np.imag(fourier)
Omega_0 = 1
delta = 3
omega = np.sqrt(Omega_0**2 + delta**2)
no_delta = 0.5*(1 - np.cos(Omega_0*t))
with_delta = 0.5 * ((Omega_0**2) / (omega))*(1 - np.cos(omega*t))

fig, axes = plt.subplots(2)
axes[0].set_title("Excited Level Population without Detuning")
axes[0].plot(t, no_delta)
axes[0].set_ylim([-0.05, 1.05])
axes[0].legend([r"$P_e$"])
# axes[0].plot(t, I)
# axes[0].legend.fancybox = True
# axes[0].set_facecolor('f0f0f0')
# -- Final --
axes[1].set_title("Excited Level Population with Detuning")
axes[1].plot(t,  with_delta)
axes[1].set_ylim([-0.05, 1.05])
axes[1].legend([r"$P_e$"])

plt.tight_layout()
plt.show()
