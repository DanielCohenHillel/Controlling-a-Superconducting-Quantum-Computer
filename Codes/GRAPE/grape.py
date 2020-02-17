# Imports
import time
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spopt
from matplotlib.pyplot import subplots
import qutip as qt
import warnings
import scipy.linalg
# warnings.filterwarnings('ignore')
a = 0.4*0
state_2 = np.array([0, 0, 1])
int_pen = 1*0


class GrapePulse:
    """
    A class used to represent a pulse to be optimized by GRAPE

    ...

    Attributes
    ----------
    :param ndarray psi_initial:
        The initial state of the qubit
    :param ndarray psi_target:
        The desired state of the qubit
    :param float total_time:
        The total time of the pulse
    :param float Ns:
        Number of time divisions, the time vector is given by t = numpy.linspace(0, total_time, Ns)
    :param ndarray base_hamiltonian:
        The hamiltonian of the system  that is not time dependent(not controled by the drives)
    :param list(ndarray) drive_hamiltonians:
        List of the time dependent hamiltonians(the list itself is not time dependent), The total, time dependent
        hamiltonian is given by
            H (t) = base_hamiltonian + drive_hamiltonian[0] * pulse0(t) + drive_hamiltonian[1] * pulse1(t) + . . .
    :param numpy.ndarray initial_pulse:
        The initial GRAPE pulse guess, the pulse matrix is of the shame [n,m] where n is the number of drive
        hamiltonians( len(drive_hamiltonians) i.e. different pulse for each hamiltonian) and m is the number of time
        steps(an attribute)
    :param boolean constraints:
        Determine weather or not to constraint the GRAPE with amplitude/slope/etc penelties
    :param float max_amp:
        The absolute maximum amplitude of the pulse, this is the hard cut-off
    :param float epsilon_soft_max:
        The soft limit on the amplitude - NOT CURRENTLY USED
    :param float lambda_band_lin:
        The strength of the slope(derivative) penalty, the bigger this is the smoother the final pulse would be
    :param boolean fix_amp_max:
        If false, the initial guess would be changed a bit by the process of the hard cut-off. Making this true will fix
        the problem but if the initial guess amplitude exceeds the maximum amplitude an error will be raised
    :param boolean print_fidelity:
        Weather or not to print the fidelity of each iteration
    """

    def __init__(self, psi_initial, psi_target,  total_time, Ns, base_hamiltonian, drive_hamiltonians,
                 initial_pulse, constraints=True, max_amp=1, epsilon_soft_max=1, lambda_band_lin=0.1,
                 lambda_amp_lin=0.03, fix_amp_max=True, print_fidelity=False, drag=False):
        try:
            # TODO: Change to auto initialize
            self.psi_initial = np.array(psi_initial)
            self.psi_target = np.array(psi_target)
            self.total_time = np.float(total_time)
            self.Ns = int(Ns)
            self.base_hamiltonian = np.array(base_hamiltonian)
            self.Nd = len(drive_hamiltonians)
            self.drive_hamiltonians = []
            for i in range(self.Nd):
                self.drive_hamiltonians.append(np.array(drive_hamiltonians[i]))
            self.initial_pulse = np.array(initial_pulse)
            self.constraints = bool(constraints)
            self.max_amp = np.float(max_amp)
            self.epsilon_soft_max = np.float(epsilon_soft_max)
            self.lambda_band_lin = np.float(lambda_band_lin)
            self.lambda_amp_lin = np.float(lambda_amp_lin)
            self.fix_amp_max = bool(fix_amp_max)
            self.print_fidelity = bool(print_fidelity)
            self.drag = bool(drag)
        except ValueError:
            raise ValueError("Invalid Input")
        except TypeError:
            raise TypeError('Invalid Input')

        self.times = np.linspace(0.0, total_time, Ns)
        self.dt = total_time / Ns
        self.dims = len(self.psi_initial)

        self._check_input()

    def _check_input(self):
        """
        Check that all the input of the GRAPE pulse are actually valid and in the right type/dimension
        :return: None
        """
        if len(self.psi_initial.shape) == 1:
            warnings.warn("Initial state of shape" + str(self.psi_initial.shape) + " is a row vector(shape == (N,))."
                          " The initial state should be a column vector(quantum state is a \"ket |>\"). \n"
                          "This will be fixed automatically but it's a bad practice to live it like so")
            self.psi_initial = np.reshape(
                self.psi_initial, (self.psi_initial.shape[0], 1))

        if len(self.psi_target.shape) == 1:
            warnings.warn("Target state of shape" + str(self.psi_target.shape) + " is a row vector(shape == (N,))."
                          " The initial state should be a column vector(quantum state is a \"ket |>\"). \n"
                          "This will be fixed automatically but it's a bad practice to live it like so")
            self.psi_target = np.reshape(
                self.psi_target, (self.psi_target.shape[0], 1))

        # Check that initial and target qubit states are actually qubit states in terms of dimensions
        init_shape = self.psi_initial.shape
        targ_shape = self.psi_target.shape
        # if not (init_shape == (2,) and targ_shape == (2,)):
        # raise ValueError('Non-pure quantum states are not currently supported in this GRAPE implementation, Make'
        # 'sure psi_initial and psi_target are of shape (2,). \npsi_initial=' + str(init_shape)
        # + '\npsi_target=' + str(targ_shape))

        # Checks if initial pulse is of correct dimensions based on the amount of drives
        if not(self.initial_pulse.flatten().shape == (self.Nd * self.Ns),):
            raise ValueError('The initial pulse must be of shape \n*amount of drive hamiltonians* by '
                             '*number of time steps*, or \n*amount of drive hamiltonians* times '
                             '*number of time steps* by 1 \n'
                             'initial_pulse is ' +
                             str(self.initial_pulse.shape) +
                             ' but it needs to be '
                             + str((self.Nd, self.Ns)))

        # Checks for valid input pulse(pulse is valid if it does not exceed the max amp) if fix_amp_max is True
        if self.fix_amp_max:
            if np.abs(np.max(self.initial_pulse)) > self.max_amp:
                raise ValueError('If \'fix_amp_max\' option is True, the initial guess MUST no exceed the maximum '
                                 'amplitude(max_amp)')

        return None

    def optimize(self):
        """
        Optimize the pulse with the GRAPE algorithm
        :return: A tuple, the first is the final pulse(same dimensions as the initial pulse) and the second
        is the final fidelity achieved by GRAPE, success is roughly if fidelity > 0.999
        """
        # Using the L-BFGS-B optimization algorithm to find the minimum of the cost function
        if self.fix_amp_max:
            result = spopt.fmin_l_bfgs_b(self.cost,
                                         np.arctanh(
                                             self.initial_pulse / self.max_amp), self.cost_gradient,
                                         factr=1e12)
        else:
            result = spopt.fmin_l_bfgs_b(self.cost, self.initial_pulse,
                                         self.cost_gradient, factr=1e9)
        result = (result[0].reshape(self.Nd,
                                    self.Ns), result[1])
        self.run_operator(result[0], show_prob=True)
        return result[0:2]

    def cost(self, in_pulse):
        """
        The cost function that the optimization algorithm seek to minimize
        :param in_pulse:
        :return: The cost function -fidelity + penalties
        """
        pulse = in_pulse.reshape(self.Nd, self.Ns)

        fid = self.run_operator(pulse)
        # print("===> ", self.constraints * self.constraint(pulse))
        fid_total = fid - self.constraints * self.constraint(pulse)
        return -fid_total

    def cost_gradient(self, in_pulse, debug_fidelity=False):
        """
        The gradient of the cost function used to minimize the cost function with the optimization algorithm
        :param in_pulse:
        :param bool debug_fidelity: Displays graphs comparing analytical and numerical fidelities for debugging
        :return: The total gradient of the cost function
        """
        # Reshape pulse the a more comfertable format [drive, time step]
        pulse = in_pulse.reshape(self.Nd, self.Ns)
        # U = self.U
        U = self.eigy_expm((1j * self.dt) * self.H(pulse))

        # -- psi_fwd --
        psi_fwd = []
        for k in range(self.Ns+1):
            if k == 0:
                psi_fwd.append(self.psi_initial)
            else:
                # print(k, "->\n", psi_fwd[-1])
                psi_fwd.append(U[k - 1] @ psi_fwd[-1])
        psi_fwd = np.array(psi_fwd)
        # print(psi_fwd.shape)
        # self.psi_fwd = psi_fwd

        # -- psi_bwd --
        # psi_bwd = np.array(
        # [np.identity(self.dims)]*(self.Ns+1), dtype=complex)
        # psi_bwd[-2] = U[-1]
        psi_bwd = np.array([np.identity(len(self.psi_initial))]
                           * (self.Ns+1), dtype=complex)

        for k in reversed(range(self.Ns)):
            psi_bwd[k] = psi_bwd[k+1]@U[k]
        phi_bwd = np.zeros([self.Ns+1, 1, len(self.psi_target)], dtype=complex)

        # Multiply the U product by <psi_target| from the left
        # for k in range(self.Ns + 1):
        #     psi_bwd[k] = (self.psi_target.conj().T @ psi_bwd[k])
        for k in range(self.Ns + 1):
            phi_bwd[k] = (self.psi_target.conj().T @ psi_bwd[k])[0]
        # print("bwd - >\n", psi_bwd[0])
        # print("f_bwd - >\n", phi_bwd[0])
        # print(k, "->\n", phi_bwd[k])

        # print("fwd-0: ", psi_fwd[0])
        # print("bwd-N: ", phi_bwd[-1])
        dc = np.ndarray(self.Ns * self.Nd, dtype=complex)
        for i, H_k in enumerate(self.drive_hamiltonians):
            for k in range(self.Ns):
                # print((phi_bwd[k+1] @ H_k @ psi_fwd[k+1])[0, 0].shape)
                dc[k + i*self.Ns] = (phi_bwd[k+1] @ H_k @ psi_fwd[k+1])[0]

        prod = np.identity(len(self.psi_initial))
        for k in range(self.Ns):
            prod = U[k] @ prod
        c = self.psi_target.conj().T @ prod @ self.psi_initial

        # --- Gradient ---
        # Before constraints
        gradient = 2 * np.real(c * np.conjugate(1j * self.dt * dc))
        # Add constraints
        gradient -= self.constraints * \
            self.constraint_gradient(pulse).flatten()

        # --- Debugging Tool ---
        if debug_fidelity:
            c_fin_transpose = np.transpose(-gradient)
            fig, axes = subplots(3)
            axes[0].set_title("gradient")
            axes[0].plot(
                self.times, -c_fin_transpose[0:self.Ns])
            axes[0].plot(
                self.times, -c_fin_transpose[self.Ns:2 * self.Ns])

            axes[1].set_title("rough estimation of gradient")
            grad = -spopt.approx_fprime(np.ndarray.flatten(
                in_pulse), self.cost, 1e-9)
            axes[1].plot(self.times, grad[0:self.Ns])
            axes[1].plot(self.times, grad[self.Ns:])

            axes[2].set_title("QI")
            axes[2].plot(self.times, pulse[0])
            axes[2].set_title("QQ")
            axes[2].plot(self.times, pulse[1])
            print("All cose: ", np.allclose(
                grad[1: -1], c_fin_transpose[1:-1], rtol=1),)

        return -gradient

    def constraint(self, pulse):
        """
        Calculate the Lagrange multipliers to create penelties on the cost function to make 'soft' limits
        :param pulse:
        :return: The total Lagrange multipliers
        """
        constraint_total = 0
        pulse = pulse.reshape(self.Nd, self.Ns)

        # --- Amplitude ---
        amp_const = np.average(pulse**2)
        amp_const *= self.lambda_amp_lin

        # --- Bandwidth ---
        slope = pulse[:, 1:] - pulse[:, :-1]
        band_const = np.sum(slope**2)
        band_const *= self.lambda_band_lin

        # --- DRAG ---
        if self.drag:
            # -- psi_fwd --
            psi_fwd = [self.psi_initial]
            for k in range(self.Ns):
                psi_fwd.append(self.U[k - 1] @ psi_fwd[-1])

            # -- Forbidden --
            forb_const = np.sum(np.abs(state_2 @ psi_fwd)**2)
            forb_const *= a
            constraint_total = forb_const

        int_const = 0
        for i in range(self.Nd):
            int_const += np.average(pulse[i])**2
        int_const /= self.Nd

        # --- Total ---
        constraint_total += amp_const + band_const + int_const*int_pen
        # print(constraint_total)
        return constraint_total.flatten()

    def constraint_gradient(self, pulse):
        """
        Calculate the gradient of the Lagrange multipliers for the optimization algorithm
        :param pulse:
        :return: The gradient of the Lagrange multipliers
        """
        constraint_total = np.zeros(pulse.shape)
        # --- Amplitude ---
        g_amp_lin = 2 * self.lambda_amp_lin * pulse

        # --- Bandwidth ---
        g_band_lin = np.zeros(pulse.shape)

        g_band_lin[:, 0] = -2 * (pulse[:, 1] - pulse[:, 0])
        g_band_lin[:, -1] = 2*(pulse[:, -1] - pulse[:, -2])

        g_band_lin[:, 1:-1] = 4 * pulse[:, 1:-1] - \
            2 * (pulse[:, 2:] + pulse[:, :-2])

        g_band_lin *= self.lambda_band_lin

        # --- DRAG ---
        g_drag = np.zeros(pulse.shape)
        if self.drag:
            # U = self.U
            U = self.eigy_expm((1j * self.dt) * self.H(pulse))
            psi_fwd = [self.psi_initial]
            for k in range(self.Ns):
                psi_fwd.append(U[k] @ psi_fwd[-1])
            psi_fwd = np.array(psi_fwd)
            # psi_fwd = self.psi_fwd
            # -- psi_bwd --
            psi_bwd = np.array(
                [np.identity(self.dims)]*(self.Ns+1), dtype=complex)
            psi_bwd[-2] = U[-1]

            for k in reversed(range(self.Ns-1)):
                psi_bwd[k] = psi_bwd[k+1]@U[k+1]

            # -- psi_bwd^-1 --
            psi_bwd_inv = np.linalg.inv(psi_bwd)

            # -- Forbiden Constraint Gradient ---
            for k in range(self.Ns):
                for i in range(k, self.Ns):
                    for j in range(self.Nd):
                        phi = state_2 @ psi_bwd_inv[i] @ psi_bwd[k]
                        overlap = 1j * self.dt * \
                            phi @ self.drive_hamiltonians[j] @ psi_fwd[k+1]
                        cc = state_2 @ psi_fwd[i+1]

                        g_drag[j, k] += a*2 * \
                            np.real(cc * np.conjugate(overlap[0]))

        g_int = np.zeros(pulse.shape)
        for i in range(self.Nd):
            g_int += 2*np.average(pulse[i])/float(self.Ns)
        g_int /= self.Nd
        # --- Total ---
        constraint_total = g_band_lin + g_amp_lin + g_drag*a + g_int*int_pen
        return constraint_total.flatten()

    def run_operator(self, pulse, show_bloch=False, calc_fidelity=True, show_prob=False):
        """
        Runs the pulse to get the final state and fidelity after driving the qubit
        :param pulse: Drive pulse to solve for the schrodinger equation
        :param boolean show_bloch: Weather or not to show the path of the qubit along the bloch sphere
        :param boolean calc_fidelity: Weather or not to calculate and return the fidelity
        :param boolean show_prob: Weather or not to show the probabilities of each state over the pulse
        :return:
        """
        prod = np.identity(self.dims)
        U = self.eigy_expm((1j * self.dt) * self.H(pulse))
        self.U = U
        if show_bloch:
            B = qt.Bloch()
            B.add_states(qt.Qobj(self.psi_initial))
            B.add_states(qt.Qobj(self.psi_target))

        if show_prob:
            Ui = self.eigy_expm((1j * self.dt) *
                                self.H(np.reshape(self.initial_pulse, (self.Nd, self.Ns))))
            initial_prob = np.zeros(
                [self.dims, self.Ns])
            final_prob = np.zeros(
                [self.dims, self.Ns])
            fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
            ax1.set_title("Initial")
            ax2.set_title("Final")
            leg = []
            wig = []
            for i in range(self.dims):
                prod = np.identity(self.dims)
                prodi = np.identity(self.dims)
                leg.append("|" + str(i) + ">")
                for k in range(self.Ns):
                    prodi = Ui[k] @ prodi
                    prod = U[k] @ prod
                    initial_prob[i, k] = np.abs(
                        (prodi @ self.psi_initial)[i]) ** 2
                    final_prob[i, k] = np.abs(
                        (prod @ self.psi_initial)[i]) ** 2
                ax1.plot(self.times, initial_prob[i, :])
                ax2.plot(self.times, final_prob[i, :])
            ax1.legend(leg)
            psi_final = prod @ self.psi_initial
            return psi_final
        else:
            for k in range(self.Ns):
                prod = U[k] @ prod
                if show_bloch:
                    print((prod @ self.psi_initial).shape)
                    B.add_points(prod @ self.psi_initial)
                    time.sleep(0.1)
            psi_final = prod @ self.psi_initial
        if show_bloch:
            B.show()
        if calc_fidelity:
            c = (self.psi_target.conj().T @ psi_final)[0, 0]
            self.c = c
            fid = np.abs(c) ** 2
            print("\n-> Fidelity: ", fid) if self.print_fidelity else None
            return fid
        return psi_final

    def H(self, pulse):
        """
        Get the total hamiltonian(base + drive hamiltonians) of a given pulse at a given index
        :param pulse: Pulse to calculate the drive hamiltonian
        :param k: Index of the hamiltonian (from 0 to Ns)
        :return: The total hamiltonian at index k
        """
        H = np.array([self.base_hamiltonian] * self.Ns)

        for i, H_k in enumerate(self.drive_hamiltonians):
            H += pulse[i, :].reshape(self.Ns, 1, 1)*H_k
        return H

    def eigy_expm(self, A, method="eigen"):
        if method == "eigen":
            vals, vects = np.linalg.eig(A)
            return np.einsum('...ik, ...k, ...kj -> ...ij',
                             vects, np.exp(vals), np.linalg.inv(vects))
        if method == "direct":
            # print("\n\n\n -> ", A.shape, "\n\n")
            for i in range(len(A)):
                A[i] = scipy.linalg.expm(A[i])
            return A
