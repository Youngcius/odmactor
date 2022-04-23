import pyvisa
from pymeasure import instruments


class LockInAmplifier(instruments.srs.SR830):
    def __init__(self, **kwargs):
        rm = pyvisa.ResourceManager()
        rs_list = rm.list_resources()
        found = False
        for rs in rs_list:
            if rs.startswith('GPIB0::') and rs.endswith('INSTR'):
                super(LockInAmplifier, self).__init__(rs, **kwargs)
                print('Lock-in Amplifier is found on {}:'.format(rs), self.id)
                found = True
                break
        if not found:
            print('Lock-in Amplifier is not found!')
