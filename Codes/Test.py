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


def H(H0, Hq_I, Hq_Q, QI, QQ, k):
    return H0 + Hq_I*QI[k] + Hq_Q*QQ[k]


print("\n---------------------------------------------------------------------\n")

h_bar = 1

# Basic operators
q = destroy(2)
qd = q.dag()

psi_initial = basis(2, 0)
psi_target = basis(2, 1)

Ns = 50
T = 5
dt = T/Ns
times = np.linspace(0.0, T, Ns)

epsilon_max = 0.1
epsilon_soft_max = 1


H0 = 1.5*qd*q
Hq_I = q + qd
Hq_Q = 1j*(q-qd)

H_d = H0
H_c = [Hq_I, Hq_Q]

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

NN = 200

result = qt.control.grape_unitary(U_targ, H_d, H_c, NN, times)

for i in range(NN):
    plt.plot(times, result.u[i,0,:])


QI = result.u[-1, 0, :]
QQ = result.u[-1, 1, :]


prod = 1
for k in range(Ns):
    U_k = (-(1j*dt)*(H0 + Hq_I*QI[k] + Hq_Q*QQ[k])).expm()
    prod = U_k*prod
psi_final = prod*psi_initial
print(ptrace(psi_target.dag()*psi_final, 0))

print(fidelity(psi_final, psi_target))
