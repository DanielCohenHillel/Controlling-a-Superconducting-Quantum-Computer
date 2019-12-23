from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qmconfig import config, HOST_IP
from time import sleep

class qmshell:
    def __init__(self):
        self.qmManager = QuantumMachinesManager(host=HOST_IP)
        self.config = config()

    # todo -> generate a new pulse operation and adding it into the config dict.
    def add_operation(self):
        pass



    def _execute(self, _func, wait_sec=5 ):
        with program() as program_object:
            I , Q , a  = ( declare(fixed) for _ in range(3) )
            _func(I, Q, a)

        job = self.qmManager.open_qm(
            self.config).execute(program_object)

        sleep(wait_sec)
        return job.get_results()

    def reset(self):

        def _rest(I, Q, A):
            pass

        self._execute( _rest, wait_sec=1)


    def __call__(self, _func, wait_sec=5):
        self.reset()
        return self._execute(_func, wait_sec)


