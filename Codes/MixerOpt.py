# Imports
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
import numpy as np
import qcodes.instrument_drivers.signal_hound.USB_SA124B
from scipy.optimize import fmin, brute
from sys import exit


def optimize():
    """
    This function optimizes the DC offsets of the quantum controller to minimize the leakage at 7GHz
    :return: Prints the ideal values of the DC offsets and correction matrix variables
    """

    # Connects to the quantum machine through the network
    qmManager = QuantumMachinesManager(host='132.77.48.243')

    # Initial guess for the best offsets
    offsets = [-0.08712208, 0.02583852]
    DC_I = offsets[0]
    DC_Q = offsets[1]

    # Initial guess for correction variables
    corvars = [0, 1]
    correction = calc_cmat(corvars[0], corvars[1])

    # Searching parameters, range of parameters for brute force, num of step in the brute force and max iteration fmin
    # OFFSETS
    nstepbruteoffset = 5  # Num of steps in the initial brute force stage(N^2)
    rangebruteoffset = [(-0.1, 0.1), (-0.1, 0.1)]  # Range to look in the inital brute force stage
    maxiterfminoffset = 10  # maximum number of iteration in the fmin stage

    # CORRECTION VARIABLES
    nstepbrutecorvars = 5  # Num of steps in the initial brute force stage(N^2)
    rangebrutecorvars = [(-0.3, 0.3), (0.6, 1.3)]  # Range to look in the inital brute force stage
    maxiterfmincorvars = 10  # maximum number of iteration in the fmin stage

    # Try to initialize the spectrum analyzer
    try:
        # Initializes the instrument
        inst = setinstrument(7e9, 25e5)

    except Exception as e:
        print("An error has occurred trying to initialize the instrument.\nThe Error:\n", e)
        exit()

    try:
        # The program that will run on the quantum machine
        with program() as prog:
            with infinite_loop_():
                play('pulse1', 'qe1')

        # The configuration of the quantum program, this is a dictionary, see documentation
        config = setconf(DC_I, DC_Q, correction)

        # Open quantum machine from configuration and force execute it
        qm1 = qmManager.open_qm(config)
        job = qm1.execute(prog, forceExecution=True)

        offsets = brute(power, rangebruteoffset, args=(corvars, qm1, inst, 'offset'), Ns= nstepbruteoffset, finish=None)

        print('\n    Initial guess for the offsets [DC_I, DC_Q]: ', offsets, "\n\n")
        # Using fmin function to find the best offsets to minimize the leakage
        xopt = fmin(power, offsets, (corvars, qm1, inst, 'offset'), maxiter=maxiterfminoffset)

    # If there's an error trying to use the "power" function
    except Exception as e:
        print("An error has occurred in the 'power' function. \nThe Error:\n", e)
        inst.close()  # Closes the spectrum analyzer
        exit()

    # Redefine the offsets to the values we found
    offsets = xopt
    print("\nOptimal offsets [DC_I, DC_Q]: " + str(offsets) + "\n\n")

    # Define the spectrum analyzer frequency to the left spike frequency
    inst.frequency(7e9 - 25e6)
    try:
        corvars = brute(powerdiffargs, rangebrutecorvars, args=(offsets, qm1, inst), Ns=nstepbrutecorvars,
                        finish=None)
        correction = calc_cmat(corvars[0], corvars[1])
        print("\n    Initial guess from brute force [th, k]: " + str(corvars) + "\n\n")

        xopt = fmin(powerdiffargs, corvars, (offsets, qm1, inst), maxiter=maxiterfmincorvars)

    except Exception as e:
        print("An error has occurred trying to find optimal correction matrix.\nThe Error:\n", e)
        exit()

    corvars = xopt
    print("\nOptimal correction variables [theta, k]: " + str(corvars))

    inst.close()  # Closes the instrument

    print("\n--------------------------------------------------------------------------------\n\n"
          "Final results: \n"
          "Offsets [DC_I, DC_Q]: " + str(offsets) +
          "\nCorrection variables [theta, k]: " + str(corvars) +
          "\n\n--------------------------------------------------------------------------------\n")


def power(offsets, corvars, qm, inst, searchfor):
    """
    This functions returns the amplitude on a specific frequency given offsets
    :param offsets: DC offsets
    :param corvars: Correction variables of the matrix, (Theta, k)
    :param qm: Quantum Manager of the program
    :param inst: The signal hound usb-SA124B
    :param searchfor: String that tells the function what is being optimized, 'offset' or 'correction'
    :return: The power of the frequency at the instrument's center frequency
    """

    # DC offsets to separate variables
    DC_I = offsets[0]
    DC_Q = offsets[1]

    # Set correction matrix
    th = corvars[0]
    k = corvars[1]
    correction = calc_cmat(th,k)
    # correction = calc_cmat(sigmoid(th)*np.pi, sigmoid(k)*0.6 + 0.7)

    # The "heart" of the function, sends the quantum machine a command to change the offsets and corrections
    if searchfor == 'offset':
        qm.set_dc_offset_by_qe('qe1', 'I', float(DC_I))
        qm.set_dc_offset_by_qe('qe1', 'Q', float(DC_Q))
    else:
        qm.set_correction('qe1', correction)

    # Read the power of the leakage frequency, this is what we want to minimize
    p = inst.get('power')

    # Prints the result of each iteration
    print("Offsets [DC_I, DC_Q]: " + str(offsets) + "  |  Correction variables [theta, k]: " + str(corvars) + " => "
          + str(p))

    return p  # Return the power of the leakage frequency


# Returns the power function but the variables are in different order for the fmin function, I will change this.
def powerdiffargs(corvars, offsets, qm, inst):
    return power(offsets, corvars, qm, inst, 'correction')


def setconf(DC_I, DC_Q, correction):
    """
    Sets the config for the quantum controller gor given DC offsets
    :param DC_I: DC offset I
    :param DC_Q: DC offset Q
    :param correction: Correction matrix
    :return: The config for the quantum machine
    """
    # Sets ports
    port_I = 3
    port_Q = 4

    LO_freq = 0
    f0 = LO_freq + 25e6

    config = {
        'version': 1,
        'controllers': {
            'con1': {
                'type': 'opx1',
                'analog': {
                    port_I: {'offset': DC_I},
                    port_Q: {'offset': DC_Q}
                }
            }
        },

        'elements': {
            'qe1': {
                'mixInputs': {
                    'I': ('con1', port_I),
                    'Q': ('con1', port_Q),
                    'mixer': 'my_mixer',
                    'lo_frequency': LO_freq
                },
                'frequency': f0,
                'operations': {
                    'pulse1': 'pulse1_in',
                }
            }
        },

        'pulses': {
            'pulse1_in': {
                'operation': 'control',
                'length': 24,
                'waveforms': {
                    'I': 'wf1',
                    'Q': 'wf2'
                }
            }

        },

        'waveforms': {
            'wf1': {
                'type': 'constant',
                'sample': 0.2
            },
            'wf2': {
                'type': 'constant',
                'sample': 0.2
            }

        },
        "mixers": {
            "my_mixer": [
                {"freq": f0, "lo_freq": LO_freq, "correction": correction}
            ]
        }
    }

    return config


def setinstrument(freq, span):
    """
    Sets up the USB-SA 124B Signal Hound spectrum analyzer and returns it as a qCoDeS instrument object
    :param freq: Center frequency of the instrument
    :param span: the span of the frequency scan
    :return: The instrument as a QCoDeS instrument
    """
    # Sets instrument drivers and connections
    path = 'C:\\Program Files\\Signal Hound\\Spike\\sa_api.dll'  # device's driver location
    mysa = qcodes.instrument_drivers.signal_hound.USB_SA124B.SignalHound_USB_SA124B('mysa', dll_path=path)

    print("\n------------------------------------------------------")
    print(mysa.get_idn())  # Prints instrument's details
    print("------------------------------------------------------\n")

    mysa.frequency(freq)  # Center of scanned region
    mysa.span(span)  # Width of scan region
    mysa.configure()

    return mysa


def calc_cmat(th, k):
    """
    Calculates the correction matrix from given angle and scalar
    :param th: The angle of the correction matrix
    :param k: The stretch of the correction matrix
    :return: The correction matrix
    """

    R = [[np.sin(th), np.cos(th)], [np.cos(th), np.sin(th)]]
    c = [[k, 0], [0, 1 / k]]
    M = np.dot(c, R).flatten().tolist()
    return M


if __name__ == '__main__':
    optimize()
