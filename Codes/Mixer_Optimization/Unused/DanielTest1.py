# Imports
import qcodes as qc
import qcodes.instrument_drivers.signal_hound.USB_SA124B
from qcodes.dataset.measurements import Measurement
from qcodes.dataset.plotting import plot_by_id


def setinstrument(freq, span):
    mysa = qcodes.instrument_drivers.signal_hound.USB_SA124B.SignalHound_USB_SA124B('mysa', dll_path='C:\\Program Files\\Signal Hound\\Spike\\sa_api.dll')  # Set instrument

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


def power(inst):
    p = inst._get_power_at_freq()
    return p


if __name__ == '__main__':
    inst = setinstrument(7e9, 1.5e7)  # Initializes the instrument

    data = scanfreqamp(inst, 2, False)  # Scan given the instruments parameters
    print(power(inst))  # Prints the power at a specified frequency
    inst.close()  # Closes the instrument at the end so we can run the code again with the same instrument name
