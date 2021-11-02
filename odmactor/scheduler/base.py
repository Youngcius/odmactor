import abc
import datetime
import time
import json
import os
import warnings
import threading
import pickle
import numpy as np
import TimeTagger as tt
from RsInstrument import RsInstrument
from odmactor.utils.asg import ASG
from typing import List, Any
import scipy.constants as C
from odmactor.utils.sequence import seq_to_str, seq_to_fig
from odmactor.utils import cal_contrast

"""
Scheduler abstract base class
"""


class Scheduler(abc.ABC):
    """
    ODMR manipulation scheduler base class
    """

    def __init__(self, *args, **kwargs):
        self._cache: Any = None
        self._data = []
        self.name = 'Base Scheduler'
        # pi pulse, for spin manipulation
        self.pi_pulse = {'freq': None, 'power': None, 'time': None}  # unit: Hz, dBm, s
        self._result = []
        self._result_detail = {}
        self._freqs = []
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
        self.time_pad = 0.0
        self.time_total = 0.0  # total time for scanning frequencies (estimated)
        self.output_dir = '../output/'
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        if 'mw_ttl' in kwargs.keys():
            self.mw_ttl = kwargs['mw_ttl']  # 1: high-level effective; 0: low level effective
        else:
            self.mw_ttl = 1  # default: high-level effective

    def set_mw_scan_freq_center_span(self, center: float, span: float, step: float, dwell: float = None):
        """
        Frequency-scanning mode with center and span frequencies to set.
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
        if self.mw_dwell != 0:
            # has set ODDR sequences
            self.time_total = self.mw_dwell * n_freqs

    def set_mw_scan_freq_start_stop(self, start: float, end: float, step: float, dwell: float = None):
        """
        Frequency-scanning mode with start and end frequencies to set.
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
        if self.mw_dwell != 0:
            self.time_total = self.mw_dwell * n_freqs

    def set_mw_freqs(self, start, end, step):
        # unit: Hz
        n_freqs = int((end - start) / step + 1)
        self._freqs = np.linspace(start, end, n_freqs).tolist()
        if self.mw_dwell != 0:
            # has set ODMR sequences
            self.time_total = self.mw_dwell * n_freqs

    def asg_connect_and_download_data(self, asg_data: List[List[float]]):
        is_connected = self._asg.connect()  # auto stop
        if is_connected == 1:
            self._asg.download_ASG_pulse_data(asg_data, [len(row) for row in asg_data])

    def configure_tagger_counting(self, apd_channel: int = None, asg_channel: int = None, reader: str = 'counter',
                                  double_readout: bool = False):
        """
        Configure asg-channel and apd-channel for ASG. For Swabian Time Tagger, channel number range: [1, 8].
        :param apd_channel: APD channel number
        :param asg_channel: ASG channel number
        :param reader: counter of specific readout type
        :param double_readout: whether use double-pulse readout
        """

        if apd_channel is not None:
            self.tagger_input['apd'] = apd_channel
        if asg_channel is not None:
            self.tagger_input['asg'] = asg_channel
        print('Current Tagger input channels:', self.tagger_input)

        # construct & execute Measurement instance
        t_ps = int(self._asg_conf['t'] / C.pico)
        N = self._asg_conf['N']

        if reader == 'counter':
            self.counter = tt.Counter(self.tagger, channels=[self.tagger_input['apd']], binwidth=t_ps, n_values=N)
        elif reader == 'cbm':
            if double_readout:
                self.counter = tt.CountBetweenMarkers(self.tagger, self.tagger_input['apd'],
                                                      begin_channel=self.tagger_input['asg'],
                                                      end_channel=-self.tagger_input['asg'],
                                                      n_values=N * 2)
            else:
                self.counter = tt.CountBetweenMarkers(self.tagger, self.tagger_input['apd'],
                                                      begin_channel=self.tagger_input['asg'],
                                                      end_channel=self.tagger_input['asg'], n_values=N)
        else:
            raise ValueError('unsupported reader (counter) type')

    @abc.abstractmethod
    def configure_odmr_seq(self, *args, **kwargs):
        pass

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

    def _get_data(self):
        self._data.append(self.counter.getData().ravel().tolist())
        self.counter.clear()

    def run(self):
        """
        1) start device
        2) acquire data timely
        """
        self._start_device()
        self._acquire_data()
        self.stop()

    def stop(self):
        self.counter.stop()
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
        self._mw_instr.write_bool('OUTPUT:STATE', True)

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

    def laser_on_seq(self, t):
        idx_laser_channel = self.channel['laser'] - 1
        self._asg_sequences[idx_laser_channel] = [t, 0]
        # self.asg.close_device()
        self.asg_connect_and_download_data(self._asg_sequences)
        self.asg.start()

    def laser_off_seq(self):
        idx_laser_channel = self.channel['laser'] - 1
        self._asg_sequences[idx_laser_channel] = [0, 0]
        # self.asg.close_device()
        self.asg_connect_and_download_data(self._asg_sequences)
        self.asg.start()

    def _conf_time_paras(self, t, N):
        """
        Configure characteristic time parameters
        :param t: total time of one ASG sequence period, unit: ns
        :param N: repetition number of ASG sequence
        """
        self._asg_conf['t'] = t * C.nano  # unit: s
        self._asg_conf['N'] = N
        self.asg_dwell = self._asg_conf['N'] * self._asg_conf['t']  # duration without padding
        self.mw_dwell = self.asg_dwell + self.time_pad * 2
        if self._freqs:
            # has set frequencies firstly
            self.time_total = self.mw_dwell * len(self._freqs)

    def _cal_counts_result(self):
        counts = [np.mean(ls) for ls in self._data]
        self._result = [self._freqs, counts]
        self._result_detail = {
            'freqs': self._freqs,
            'counts': counts,
            'origin_data': self._data
        }

    def _cal_contrasts_result(self):
        contrasts = [cal_contrast(ls) for ls in self._data]
        self._result = [self._freqs, contrasts]
        self._result_detail = {
            'freqs': self._freqs,
            'contrast': contrasts,
            'origin_data': self._data,
        }

    def save_result(self, fname: str = None):
        if not self._result or not self._result_detail:
            raise ValueError('empty result cannot be saved')
        if fname is None:
            fname = os.path.join(self.output_dir,
                                 '{}-result-{}-{}'.format(self.name, str(datetime.date.today()),
                                                          round(time.time() / 120)))
        with open(fname + '.json', 'w') as f:
            json.dump(self._result_detail, f)
        np.savetxt(fname + '.txt', self._result)
        print('result has been saved into {}'.format(fname + '.json'))

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

    @property
    def sequences_strings(self):
        return seq_to_str(self._asg_sequences)

    @property
    def sequences_figure(self):
        return seq_to_fig(self._asg_sequences)


def dBm_to_mW(dBm):
    return 10 ** (dBm / 10)


def mW_to_dBm(mW):
    return 10 * np.log10(mW)


class SimpleScheduler(Scheduler):
    """
    Simple Scheduler: single-frequency scheduling
    """

    def configure_odmr_seq(self, *args, **kwargs):
        pass

    def _start_device(self, *args, **kwargs):
        pass

    def _acquire_data(self, *args, **kwargs):
        pass


class ODMRScheduler(Scheduler):
    def __init__(self, *args, **kwargs):
        super(ODMRScheduler, self).__init__(*args, **kwargs)
        self.name = 'Base ODMR Scheduler'

    def _start_device(self):
        # 1. run ASG firstly
        self._data.clear()
        self._asg.start()

        # 2. run MW then
        self._mw_instr.write_bool('OUTPUT:STATE', True)
        if self.mw_exec_mode == 'scan-center-span' or self.mw_exec_mode == 'scan-start-stop':
            self._mw_instr.write_str('SWE:FREQ:EXEC')  # trigger the sweep
            raise warnings.warn('Using the SWEEP mode of MW', DeprecationWarning)
            # self._mw_instr.write_str('SOUR:PULM:TRIG:MODE SING')
            # self._mw_instr.write_str('SOUR:PULM:TRIG:IMM')
        print('MW on/off status:', self._mw_instr.instrument_status_checking)

    def _scan_freqs_and_get_data(self):
        """
        Sanning frequencies & getting data of Counter
        """
        for i, freq in enumerate(self._freqs):
            print('scanning freq {:.4f} GHz'.format(freq / C.giga))
            t = threading.Thread(target=self._get_data, name='thread-{}'.format(i))
            time.sleep(self.time_pad)
            time.sleep(self.asg_dwell)  # accumulate counts
            t.start()  # begin readout
            time.sleep(self.time_pad)
            t.join()

    def run_scanning(self, mw_control='on'):
        """
        1) start device
        2) acquire data timely
        """
        mw_seq_on = self._asg_sequences[self.channel['mw'] - 1]
        if mw_control == 'off':
            # mw_seq_off = [0, sum(mw_seq_on)]
            self._asg_sequences[self.channel['mw'] - 1] = [0, 0]
            self.asg_connect_and_download_data(self._asg_sequences)
        elif mw_control == 'on':
            pass
        else:
            raise ValueError('unsupported MW control parameter (should be "or" or "off"')
        print('Begin to run {}. Frequency: {:.4f} - {:.4f} GHz.'.format(self.name, self._freqs[0], self._freqs[-1]))
        print('t: {:.2f} ns, N: {}, T: {:.2f} s, n_freqs: {}'.format(self._asg_conf['t'] / C.nano, self._asg_conf['N'],
                                                                     self.mw_dwell, len(self._freqs)))
        print('Estimated total running time: {:.2f} s'.format(self.time_total))
        self._start_device()
        self._acquire_data()
        self.stop()

        # 恢复微波的ASG的MW通道为 on
        if mw_control == 'off':
            self._asg_sequences[self.channel['mw'] - 1] = mw_seq_on
            self.asg_connect_and_download_data(self._asg_sequences)

    @abc.abstractmethod
    def configure_odmr_seq(self, *args, **kwargs):
        pass

    @abc.abstractmethod
    def _acquire_data(self, *args, **kwargs):
        """
        Acquire data in real time
        """
        raise NotImplemented


class SpinControlScheduler(Scheduler):
    """
    Spin control abstract class
    """

    def __init__(self, *args, **kwargs):
        super(SpinControlScheduler, self).__init__(*args, **kwargs)
        self.name = 'Base Spin Control Scheduler'
