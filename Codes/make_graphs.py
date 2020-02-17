import numpy as np
import qutip as qt
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')

N = 1000
W = 50
x = np.linspace(-1, 6, N)
w = 3
sig = 1
gauss = np.exp((-(x-w)**2)/0.3)

T = 11
Nt = 1000
t = np.linspace(-T/2, T/2, Nt)
fourier = np.exp((-t**2)/4 - 1j*w*t)

plt.plot(t, np.real(fourier))
plt.plot(t, np.imag(fourier))
# plt.plot(t, np.exp((-t**2)/4), '--')
R = np.real(fourier)
I = np.imag(fourier)

fig, axes = plt.subplots(1, 2)
axes[0].set_title("Pulses in Time Space")
axes[0].plot(t, R)
axes[0].plot(t, I)
# axes[0].legend.fancybox = True
axes[0].legend([r"$\epsilon_I$", r'$\epsilon_Q$'])
# axes[0].set_facecolor('f0f0f0')
# -- Final --
axes[1].set_title("Pulses in Frequency Space")
axes[1].plot(x,  gauss)

plt.tight_layout()
plt.show()
