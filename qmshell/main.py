from qmshell import qmshell
from qm.qua import *
from matplotlib import pyplot as plt

if __name__ == '__main__':

    def _func(I , Q, a):

        #play('pulse1' * amp(10), 'qe1')
        with infinite_loop_():
            # play('pulse1' * amp(10), 'qe1')
            measure('meas_pulse_in' * amp(10), 'qe1', 'samples', ('cos_weights', Q), ('sin_weights', I))
            save(I, 'I')
            save(Q, 'Q')
            # save(a, 'a')

    _shell = qmshell()
    results = _shell(_func, wait_sec=2)

    Idata = results.variable_results.I.data
    Qdata = results.variable_results.Q.data

    plt.scatter(Idata, Qdata)
    plt.show()