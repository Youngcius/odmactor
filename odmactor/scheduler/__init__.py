import abc
import datetime
import os
import pickle
import TimeTagger as tt
import numpy as np
# from scheduler.continuity import CWScheduler
# from scheduler.pulse import PulseScheduler
# from scheduler.ramsey import RamseyScheduler
from RsInstrument import RsInstrument
from odmactor.utils.asg import ASG
from typing import List, Any

"""
Priority to be implemented: CW, Pulse, Ramsey, Rabi
TODO: DD, Hahn echo, T1 relaxation
Reference: https://arxiv.org/pdf/1910.00061.pdf
"""


class Scheduler(abc.ABC):
    """
    ODMR manipulation scheduler base class
    """

    def __init__(self, output_dir: str = '../output/', *args, **kwargs):
        self._cache: Any = None
        self.name = 'ODMR Base Scheduler'
        # pi pulse, for spin manipulation
        self.pi_pulse = {'freq': None, 'power': None, 'time': None}  # unit: Hz, dBm, s
        self._result = {}
        self._result_detail = {}
        self._freqs = np.array([])
        self._asg_sequences = []
        self.reset_asg_sequence()
        self._asg_conf = {'t': None, 'N': None}
        self._configuration = {}
        self._laser_control = True
        self._mw_instr = RsInstrument('USB0::0x0AAD::0x0054::104174::INSTR', True, True)
        self._asg = ASG()
        self.mw_exec_mode = ''
        self.mw_exec_modes_optional = {'scan-center-span', 'scan-start-stop'}
        self.channel = {'laser': 1, 'mw': 2, 'apd': 3, 'tagger': 5}
        self.tagger_input = {'apd': 1, 'asg': 2}
        self.tagger = tt.createTimeTagger()
        self.counter: tt.IteratorBase = None
        # properties or method for debugging
        self.sync_delay = 0.0
        self.mw_dwell = 0.0
        self.asg_dwell = 0.0
        self.time_pad = 1.0
        self.time_total = 0.0  # total time for scanning frequencies
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        if 'mw_ttl' in kwargs.keys():
            self.mw_ttl = kwargs['mw_ttl']  # 1: high level; 0: low level
        else:
            self.mw_ttl = 1

    def set_mw_scan_freq_center_span(self, center: float, span: float, step: float, dwell: float = None):
        """
        frequency-scanning mode with center and span frequencies to set
        :param center: center frequency, unit: Hz
        :param span: span frequency, unit: Hz
        :param step: each step frequency, unit: Hz
        :param dwell: dwelling time for each scanning point, unit: s
        """
        self.mw_exec_mode = 'scan-center-span'
        self._mw_instr.write_float('FREQ:CENT', center)
        self._mw_instr.write_float('FREQ:SPAN', span)
        self._mw_instr.write_float('SOUR:SWE:FREQ:STEP:LIN', step)  # define the step
        if dwell is not None:
            self.mw_dwell = dwell
        self._mw_instr.write_float('SWE:DWEL', self.mw_dwell)
        ########################
        # 2021-10-15 modified
        # self._mw_instr.write_str('SOUR:SWE:FREQ:MODE AUTO')
        # self._mw_instr.write_str('TRIG:FSW:SOUR EXT')  # external trigger
        ########################
        self._mw_instr.write_str('TRIG:FSW:SOUR SING')

        self._mw_instr.write_str('SOUR:FREQ:MODE SWE')  # define sweep mode
        n_freqs = int(span / step + 1)
        self._freqs = np.linspace(center - span / 2, center + span / 2, n_freqs).tolist()
        self.time_total = self.mw_dwell * n_freqs

    def set_mw_scan_freq_start_stop(self, start: float, end: float, step: float, dwell: float = None):
        """
        frequency-scanning mode with start and end frequencies to set
        :param start: start frequency, unit: Hz
        :param end: end frequency, unit: Hz
        :param step: each step frequency, unit: Hz
        :param dwell: dwelling time for each scanning point, unit: s, default: 1s
        """

        self.mw_exec_mode = 'scan-start-stop'
        self._mw_instr.write_float('FREQ:STAR', start)
        self._mw_instr.write_float('FREQ:STOP', end)
        self._mw_instr.write_float('SOUR:SWE:FREQ:STEP:LIN', step)  # define the step
        if dwell is not None:
            self.mw_dwell = dwell  # unit: s
        self._mw_instr.write_float('SWE:DWEL', self.mw_dwell)
        ########################
        # 2021-10-15 modified
        # self._mw_instr.write_str('SOUR:SWE:FREQ:MODE AUTO')
        # self._mw_instr.write_str('TRIG:FSW:SOUR EXT')  # external trigger
        ########################
        # 2021-10-21 modified
        self._mw_instr.write_str('TRIG:FSW:SOUR SING')

        self._mw_instr.write_str('SOUR:FREQ:MODE SWE')  # define sweep mode
        n_freqs = int((end - start) / step + 1)
        self._freqs = np.linspace(start, end, n_freqs).tolist()
        self.time_total = self.mw_dwell * n_freqs

    def set_mw_freqs(self, start, end, step):
        # unit: Hz
        n_freqs = int((end - start) / step + 1)
        self._freqs = np.linspace(start, end, n_freqs).tolist()
        self.time_total = self.mw_dwell * n_freqs
        self._mw_instr.write_bool('OUTPUT:STATE', True)

    def asg_connect_and_download_data(self, asg_data: List[List[float]]):
        is_connected = self._asg.connect()
        if is_connected == 1:
            self._asg.download_ASG_pulse_data(asg_data, [len(row) for row in asg_data])

    def configure_tagger_counting(self, apd_channel: int = None, asg_channel: int = None):
        """
        Configure asg-channel and apd-channel for ASG. For Swabian Time Tagger, channel number range: [1, 8].
        :param apd_channel: APD channel number
        :param asg_channel: ASG channel number
        """
        if apd_channel is not None:
            self.tagger_input['apd'] = apd_channel
        if asg_channel is not None:
            self.tagger_input['asg'] = asg_channel
        print('Current Tagger input channels:', self.tagger_input)

        # execute Measurement instance
        self.counter = tt.CountBetweenMarkers(self.tagger, self.tagger_input['apd'],
                                              begin_channel=self.tagger_input['asg'],
                                              end_channel=-self.tagger_input['asg'],
                                              n_values=int(self._asg_conf['N'] * self._asg_conf['t']))

    @abc.abstractmethod
    def _start_device(self, *args, **kwargs):
        """
        Start device: MW, ASG; Execute Measurement instance.
        """
        raise NotImplemented

    @abc.abstractmethod
    def _acquire_data(self, *args, **kwargs):
        """
        Acquire data in real time
        """
        raise NotImplemented

    def run(self):
        """
        1) start device
        2) acquire data timely
        """
        self._start_device()
        self._acquire_data()
        self.stop()

    def stop(self):
        self._asg.stop()
        # self._mw_instr.write_str('STOP')
        # TODO: how to stop the microwave

        print('Scheduling process has stopped')

    def close(self):
        self._asg.close_device()
        self._mw_instr.close()
        tt.freeTimeTagger(self.tagger)
        print('All instrument resources has been released')

    def configure_mw_paras(self, power: float = None, freq: float = None, regulate_pi: bool = False, *args, **kwargs):
        # TODO: check 'float' type
        if power is not None:
            self._mw_instr.write_float('POW', power)
            if regulate_pi:  # regulate time duration based on MW power
                self._regulate_pi_pulse(power=power)  # by power
        if freq is not None:
            self._mw_instr.write_float('FREQUENCY', freq)

    def _regulate_pi_pulse(self, power: float = None, time: float = None):
        # if by not in {'power', 'time'}:
        #     raise TypeError('Please input a supported parameter: "by" should be one of {"power", "time"}')
        if not any(self.pi_pulse.values()):  # means it has no pi-pulse parameters
            raise TypeError('No pi pulse parameters to be regulated')
        else:
            power_ori = self.pi_pulse['power']  # unit: dBm
            time_ori = self.pi_pulse['time']
            if power is not None:
                # calculate new "time"
                time = np.sqrt(10 ** ((power_ori - power) / 10)) * time_ori
            elif time is not None:
                # calculate new "power"
                power = mW_to_dBm((time_ori / time) ** 2 * dBm_to_mW(power_ori))
            self.pi_pulse['power'] = power
            self.pi_pulse['time'] = time

    def reset_asg_sequence(self):
        self._asg_sequences = [[0, 0] for i in range(8)]

    def save_configuration(self, fname: str = None):
        """
        Save configuration parameters into a disk file
        :param fname: file name for saving
        """
        if fname is None:
            fname = '{} {}'.format(self.name, datetime.date.today())
        with open(fname, 'wr') as f:
            pickle.dump(self, f)
        print('object has been save to {}'.format(fname))

    def __str__(self):
        return self.name

    @property
    def laser_control(self):
        return self._laser_control

    @property
    def mw_instr(self):
        return self._mw_instr

    @mw_instr.setter
    def mw_instr(self, value: RsInstrument):
        self._mw_instr = value

    @property
    def asg(self):
        return self._asg

    @asg.setter
    def asg(self, value: ASG):
        self.asg = value

    @property
    def result(self):
        return self._result

    @property
    def result_detail(self):
        return self._result_detail


def dBm_to_mW(dBm):
    return 10 ** (dBm / 10)


def mW_to_dBm(mW):
    return 10 * np.log10(mW)
