import time
from RsInstrument import RsInstrument


class Microwave(RsInstrument):
    """
    Microwave class based on RsInstrument (vendor: R&S)
    """

    def __init__(self):
        super(Microwave, self).__init__('USB0::0x0AAD::0x0054::104174::INSTR', True, True)

    def set_frequency(self, freq):
        super(Microwave, self).write_float('FREQUENCY', freq)

    def set_power(self, power):
        super(Microwave, self).write_float('POW', power)

    def run_given_time(self, duration):
        self.start()
        time.sleep(duration)
        self.stop()

    def connect(self, force_close: bool = False) -> bool:
        return super(Microwave, self).reconnect(force_close)

    def start(self):
        self.write_bool('OUTPUT:STATE', True)

    def stop(self):
        self.write_bool('OUTPUT:STATE', False)

    def close(self):
        super(Microwave, self).close()
