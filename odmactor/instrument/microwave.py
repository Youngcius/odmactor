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

# # Simple example on how to use the RsInstrument module for remote-controlling yor VISA instrument
# # Preconditions:
# # - Installed RsInstrument Python module (see the attached RsInstrument_PythonModule folder Readme.txt)
# # - Installed VISA e.g. R&S Visa 5.12.x or newer

# from RsInstrument.RsInstrument import RsInstrument

# resource_string_1 = 'TCPIP::192.168.2.101::INSTR'  # Standard LAN connection (also called VXI-11)
# resource_string_2 = 'TCPIP::192.168.2.101::hislip0'  # Hi-Speed LAN connection - see 1MA208
# resource_string_3 = 'GPIB::20::INSTR'  # GPIB Connection
# resource_string_4 = 'USB::0x0AAD::0x0119::022019943::INSTR'  # USB-TMC (Test and Measurement Class)
# resource_string_5 = 'RSNRP::0x0095::104015::INSTR'  # R&S Powersensor NRP-Z86
# instr = RsInstrument(resource_string_1, True, False)

# idn = instr.query_str('*IDN?')
# print(f"\nHello, I am: '{idn}'")
# print(f'RsInstrument driver version: {instr.driver_version}')
# print(f'Visa manufacturer: {instr.visa_manufacturer}')
# print(f'Instrument full name: {instr.full_instrument_model_name}')
# print(f'Instrument installed options: {",".join(instr.instrument_options)}')

# # Close the session
# instr.close()
