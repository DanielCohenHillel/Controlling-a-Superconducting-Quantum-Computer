import numpy as np




# globals
HOST_IP = '132.77.48.243'
PORT_Q = 2
PORT_I = 1




# ---------------------------
# Daniel changed this line, You can chane it back to the line above
# offsets = [-0.08949547, 0.02586283]
# offsets = [-0.08232403, 0.0194166]
offsets = [-0.08302877, 0.01893423]
DC_I = offsets[0]
DC_Q = offsets[1]


LO_freq = 0
f0 = LO_freq + 25e6

#corvars = [2.16005713e-04, 1.4]
corvars = [-0.10308351, 0.89446135]
#corvars = [-0.11507042, 0.89980469]
th = corvars[0]
k = corvars[1]

R = [[np.sin(th),np.cos(th)], [np.cos(th), np.sin(th)]]
c = [[k, 0], [0 , 1/k]]
M = np.dot(c, R).flatten().tolist()
# ------------------------------


def config(port_i=PORT_I, port_q=PORT_Q):
    return {
        'version': 1,
        'controllers': {
            'con1': {
                'type': 'opx1',
                'analog_outputs': {
                    port_i: {'offset': DC_I},
                    port_q: {'offset': DC_Q}
                }
            }
        },

        'elements': {
            'qe1': {
                'mixInputs': {
                    'I': ('con1', port_i),
                    'Q': ('con1', port_q),
                    'mixer': 'my_mixer',
                    'lo_frequency': LO_freq
                },
                'intermediate_frequency': f0,
                'operations': {
                    'pulse1': 'pulse1_in',
                    'meas_pulse': 'meas_pulse_in'
                },
                # todo | just a patch
                'time_of_flight': 2000,
                'smearing': 0,
                'outputs': {
                    'out1': ('con1', 1)
                },
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
            },
            'meas_pulse_in': {
                'operation': 'measurement',
                'length': 200,
                'waveforms': {
                    'I': 'wf1',
                    'Q': 'wf2'
                },
                'integration_weights': {
                    'cos_weights': 'cos_weights',
                    'sin_weights': 'sin_weights',
                },
                'digital_marker': 'marker1'
            },
        },

        'waveforms': {
            'wf1': {
                'type': 'constant',
                'sample': 0.1
            },
            'wf2': {
                'type': 'constant',
                'sample': 0.1
            }

        },

        'digital_waveforms': {
            'marker1': {
                'samples': [(1, 4), (0, 2), (1, 1), (1, 0)]
            }
        },

        "integration_weights": {
            "cos_weights": {
                'cosine': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                'sine': [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
                       4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
                       4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
                       4.0, 4.0, 4.0, 4.0, 4.0, 4.0]
            },
            "sin_weights": {
                'cosine': [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
                           4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
                           4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
                           4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
                'sine': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            }
        },

        "mixers": {
            "my_mixer": [
                {"intermediate_frequency": f0, "lo_frequency": LO_freq, "correction": [1.0, 0.0, 0.0, 1.0]}
            ]
        }

    }


