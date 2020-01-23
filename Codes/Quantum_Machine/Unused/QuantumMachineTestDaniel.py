# Imports
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
import qm as qm
import time
import numpy as np
import qcodes as qc
import qcodes.instrument_drivers.signal_hound.USB_SA124B
from qcodes.dataset.measurements import Measurement
from qcodes.dataset.plotting import plot_by_id
from scipy.optimize import fmin
import matplotlib
import random


# This function optimizes the DC offsets of the quantum controller to minimize the leakage at 7GHz
def optimize():
    # Connects to the quantum machine through the network
    qmManager = QuantumMachinesManager(host='132.77.48.245')

    # Initial guess for the best offsets
    offsets = [-0.10059273, 0.1322426]
    DC_I = offsets[0]
    DC_Q = offsets[1]

    # Initializes the instrument
    inst = setinstrument(7e9-25e6, 250e3)

    try:
        # Using fmin function to find the best offsets to minimize the leakage
        xopt = fmin(power, 0.2, (qmManager, inst), maxiter=10)

    # If there's an error trying to use the "power" function
    finally:
        print("An error has occurred trying to find to minimum of the power function")
        inst.close()  # Closes the spectrum analyzer
    print(str(xopt) + " ==> " + str(power(xopt, qmManager, inst)))  # Print the results

    # inst.span(2e8)
    # scanfreqamp(inst, 20, True)
    inst.close()  # Closes the instrument


# This functions returns the amplitude on a specific frequency given offsets
def power(th, qmManager, inst):
    offsets = [-0.10059273, 0.1322426]
    DC_I = offsets[0]
    DC_Q = offsets[1]

    config = setconf(th)

    qm1 = qmManager.open_qm(config)
    #  print(qm1.__getattribute__('id'))
    #  print(qmManager.get_qm('').set_dc_offset_by_qe())

    # Plays pulse infinitely
    with program() as prog:
        with infinite_loop_():
            play('pulse1', 'qe1')

    job = qm1.execute(prog, forceExecution=True)

    # inst.frequency(7e9)

    p = inst._get_power_at_freq()

    '''
    inst.frequency(7e9 - 25e6)
    print("************************")
    print("==>  " + str(inst._get_power_at_freq()) + "  <==")
    print("************************")
    '''

    print(str(th) + " => " + str(p))
    qm1.close()
    job.halt()
    return p


# Sets the config for the quantum controller gor given DC offsets
def setconf(th):
    offsets = [-0.10059273, 0.1322426]
    DC_I = offsets[0]
    DC_Q = offsets[1]

    # Sets ports
    port_I = 9
    port_Q = 10

    LO_freq = 0
    f0 = LO_freq + 25e6

   # th = 0.1  # rotation angle
    R = [[np.cos(th), np.sin(th)], [np.sin(th), -np.cos(th)]]  # Rotation matrix
    c = [[1, 0], [0, 0.8]]

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
                'sample': 0.1 * 0
            },
            'wf2': {
                'type': 'constant',
                'sample': 0.05
            }

        },
        "mixers": {
            "my_mixer": [
                {"freq": f0, "lo_freq": LO_freq, "correction": [random.randrange(1, 100)/50, random.randrange(1, 100)/50]*2}
            ]
        }
    }

    return config


# Sets up the USB-SA 124B Signal Hound spectrum analyzer and returns it as a qCoDeS instrument object
def setinstrument(freq, span):
    # Sets instrument drivers and connections
    path = 'C:\\Program Files\\Signal Hound\\Spike\\sa_api.dll'  # device's driver location
    mysa = qcodes.instrument_drivers.signal_hound.USB_SA124B.SignalHound_USB_SA124B('mysa', dll_path=path)

    print("-------------------------------------")
    print(mysa.get_idn())  # Prints instrument's details
    print("-------------------------------------")

    mysa.frequency(freq)  # Center of scanned region
    mysa.span(span)  # Width of scan region

    return mysa


def scanfreqamp(inst, avg, plot = False):
    inst.avg(avg)  # Number of traces to avarage over, the "quality" of the measurement
    meas = Measurement()  # Creates new measurement
    meas.register_parameter(inst.freq_sweep)
    with meas.run() as datasaver:
        datasaver.add_result((inst.frequency_axis, inst.frequency_axis.get()),
                             (inst.freq_sweep, inst.freq_sweep.get(),))

        dataid = datasaver.run_id

    if plot:
        plot_by_id(dataid)  # Plots the measurement
    return datasaver


if __name__ == '__main__':
    optimize()



# Some results

# [ 0.08293153  0.50117324 -0.04978425  0.15800019] => 14.268719673156738

# [-0.10093183  0.13161536 -0.0953647   0.10385427] => -87.81173706054688

# [-0.10106181  0.13327702]

# [-0.10101066  0.13224149]

# [-0.1010631   0.13629986] ==> -80.63866424560547

# [-0.1006443   0.13447224] ==> -69.32542419433594

# [-0.10065907  0.13449088] => -83.85325622558594

# [-0.10062767  0.1419932 ] => -87.35037994384766

# [-0.10062767  0.1419932 ] ==> -50.07796859741211

#[-0.10066698  0.13336822] ==> -75.14099884033203

#[-0.1005527   0.13180124] => -85.38912200927734

# [-0.10059273  0.1322426 ] ==> -88.18570709228516
