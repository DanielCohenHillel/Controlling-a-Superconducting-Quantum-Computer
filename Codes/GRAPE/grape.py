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
# a = 0.4*0
state_2 = np.array([0, 0, 1])
int_pen = 1*0
n_ph = 5


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
                 lambda_amp_lin=0.03, lambda_drag=0.1, fix_amp_max=True, print_fidelity=False, drag=False):
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
            self.lambda_drag = np.float(lambda_drag)
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

        # Check if initial and target states match in dimensions
        init_shape = self.psi_initial.shape
        targ_shape = self.psi_target.shape
        if not (init_shape == targ_shape):
            raise ValueError("Initial state dimensions " + str(init_shape)
                             + " must match target state dimensions " + str(targ_shape))

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
                                         self.cost_gradient, factr=1e10)
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

        U = self.eigy_expm((1j * self.dt) * self.H(pulse))

        # -- psi_fwd --
        psi_fwd = []
        for k in range(self.Ns+1):
            if k == 0:
                psi_fwd.append(self.psi_initial)
            else:
                psi_fwd.append(U[k - 1] @ psi_fwd[-1])
        psi_fwd = np.array(psi_fwd)

        # -- psi_bwd --
        psi_bwd = np.array([np.identity(len(self.psi_initial))]
                           * (self.Ns+1), dtype=complex)

        for k in reversed(range(self.Ns)):
            psi_bwd[k] = psi_bwd[k+1]@U[k]

        # -- phi_bwd --
        phi_bwd = np.zeros([self.Ns+1, 1, len(self.psi_target)], dtype=complex)

        # Multiply psi_bwd by <psi_target| from the left into phi_bwd
        for k in range(self.Ns + 1):
            phi_bwd[k] = (self.psi_target.conj().T @ psi_bwd[k])[0]

        # -- Overlap Gradient --
        dc = np.ndarray(self.Ns * self.Nd, dtype=complex)
        for i, H_k in enumerate(self.drive_hamiltonians):
            for k in range(self.Ns):
                dc[k + i*self.Ns] = (phi_bwd[k+1] @ H_k @ psi_fwd[k+1])[0]

        # -- Overlap --
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

        itime = time.time()
        # --- DRAG ---
        if self.drag:
            # -- psi_fwd --
            psi_fwd = [self.psi_initial]
            for k in range(self.Ns):
                psi_fwd.append(self.U[k - 1] @ psi_fwd[-1])

            # -- Forbidden --
            forb_const = np.sum(np.abs(state_2 @ psi_fwd)**2)
            forb_const *= self.lambda_drag
            constraint_total = forb_const
        print("drag time: ", time.time() - itime)
        # --- Total ---
        constraint_total += amp_const + band_const
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

        itime = time.time()
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

                        g_drag[j, k] += self.lambda_drag*2 * \
                            np.real(cc * np.conjugate(overlap[0]))
            print("drag grad time ", time.time() - itime)

        # --- Total ---
        constraint_total = g_band_lin + g_amp_lin + g_drag
        return constraint_total.flatten()

    def cav_cost(self, pulse, kind='cost'):
        sum_fid = 0

        shape_cav = [[2, int(self.dims/2)], [1, 1]]

        qubit_init = np.array([[0]]*shape_cav[0][0])
        qubit_target = np.array([[0]]*shape_cav[0][0])

        cavity_init = np.array([[0]]*shape_cav[0][1])
        cavity_target = np.array([[0]]*shape_cav[0][1])

        qubit_qobj_init = qt.Qobj(self.psi_initial, shape_cav).ptrace(0)
        qubit_qobj_target = qt.Qobj(self.psi_target, shape_cav).ptrace(0)

        cavity_qobj_init = qt.Qobj(self.psi_initial, shape_cav).ptrace(1)
        cavity_qobj_target = qt.Qobj(self.psi_target, shape_cav).ptrace(1)
        for i in range(shape_cav[0][0]):
            qubit_init[i] = np.abs(qubit_qobj_init[i][0][i])
            qubit_target[i] = np.abs(qubit_qobj_target[i][0][i])
        for i in range(shape_cav[0][1]):
            cavity_init[i] = np.abs(cavity_qobj_init[i][0][i])
            cavity_target[i] = np.abs(cavity_qobj_target[i][0][i])

        og_cavity_init = cavity_init
        og_cavity_target = cavity_target

        og_base_hamiltonian = self.base_hamiltonian
        og_drive_hamiltonians = self.drive_hamiltonians

        w_a = 5.66*0.5
        w_c = 4.5
        chi = 2
        alpha = w_a*0.05
        K = w_c*10**(-3)
        chitag = chi*10**(-2)
        for i in range(n_ph):
            cavity_init = np.array([np.append(cavity_init, [0])]).T
            cavity_target = np.array([np.append(cavity_target, [0])]).T

            self.psi_initial = np.kron(qubit_init, cavity_init)
            self.psi_target = np.kron(qubit_target, cavity_target)
            self.dims = len(self.psi_initial)

            qubit_levels = len(qubit_init)
            cavity_levels = len(cavity_init)

            a = qt.tensor(qt.destroy(qubit_levels), qt.qeye(cavity_levels))
            ad = a.dag()

            c = qt.tensor(qt.qeye(qubit_levels), qt.destroy(cavity_levels))
            cd = c.dag()

            Ha = w_a*ad*a + (alpha/2)*ad*ad*a*a
            Hc = w_c*cd*c + (K/2)*cd*cd*c*c
            Hi = chi*cd*c*ad*a + chitag*cd*cd*c*c*ad*a
            H0 = np.array(Ha + Hc + Hi)

            Ha_I = np.array(a + ad)
            Ha_Q = np.array(1j*(a - ad))

            Hc_I = np.array(c + cd)
            Hc_Q = np.array(1j*(c - cd))

            drive_hamiltonians = [Ha_I, Ha_Q, Hc_I, Hc_Q]

            self.drive_hamiltonians = drive_hamiltonians
            self.base_hamiltonian = H0

            if kind == 'cost':
                sum_fid += self.cost(pulse)
            else:
                sum_fid += self.cost_gradient(pulse)
        cavity_init = og_cavity_init
        cavity_target = og_cavity_target
        self.psi_initial = np.kron(qubit_init, cavity_init)
        self.psi_target = np.kron(qubit_target, cavity_target)
        self.dims = len(self.psi_initial)

        self.base_hamiltonian = og_base_hamiltonian
        self.drive_hamiltonians = og_drive_hamiltonians

        return sum_fid

    def cost_gradient_cav(self, pulse):
        return self.cav_cost(pulse, kind='grad')

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
            ax1.set_title("Initial Level Population Over Time")
            ax2.set_title("Final Level PopulationTime")

            ax1.set_xlabel("Time")
            ax1.set_ylabel("Amplitude")
            ax2.set_xlabel("Time")
            ax2.set_ylabel("Amplitude")
            leg = []
            wig = []
            N_cav = int(self.dims/2)
            for i in range(self.dims):
                prod = np.identity(self.dims)
                prodi = np.identity(self.dims)
                if i < N_cav:
                    leg.append(r'$|g> \otimes |{}>$'.format(str(i)))
                else:
                    leg.append(r'$|e> \otimes |{}>$'.format(str(i-N_cav)))

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
            states_qobjs = []
            for k in range(self.Ns):
                prod = U[k] @ prod
                if show_bloch:
                    states_qobjs.append(qt.Qobj(prod @ self.psi_initial))
            psi_final = prod @ self.psi_initial
        if show_bloch:
            B.add_states(states_qobjs, kind='point')
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
            for i in range(len(A)):
                A[i] = scipy.linalg.expm(A[i])
            return A
