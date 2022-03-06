"""
Scheduler abstract base class
"""
import abc
import datetime
import time
import uuid
import json
import os
import pickle
import threading
import numpy as np
import scipy.constants as C
import TimeTagger as tt
from RsInstrument import RsInstrument
from odmactor.utils.asg import ASG
from odmactor.utils import dBm_to_mW, mW_to_dBm
from typing import List, Any, Optional
from odmactor.utils.sequence import seq_to_str, seq_to_fig
from matplotlib.figure import Figure


class Scheduler(abc.ABC):
    """
    ODMR manipulation scheduler base class
    """

    def __init__(self, *args, **kwargs):
        self._cache: Any = None
        self._data = []
        self._data_ref = []
        self.name = 'Base Scheduler'
        # pi pulse, for spin manipulation
        self.pi_pulse = {'freq': None, 'power': None, 'time': None}  # unit: Hz, dBm, s
        self._result = []
        self._result_detail = {}
        self._freqs = []
        self._times = []
        self._asg_sequences = []
        self.reset_asg_sequence()
        self._asg_conf = {'t': 0, 'N': 0}  # to calculate self.asg_dwell = N * t
        self._mw_conf = {'freq': C.giga, 'power': 0}  # current MW parameter settings
        self._configuration = {}
        self._laser_control = True
        self.two_pulse_readout = False  # whether use double-pulse readout
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

        if 'with_ref' in kwargs.keys():
            self.with_ref = kwargs['with_ref']
        else:
            self.with_ref = False
        if 'epoch_omit' in kwargs.keys():
            self.epoch_omit = kwargs['epoch_omit']
        else:
            self.epoch_omit = 0

    def asg_connect_and_download_data(self, asg_data: List[List[float]]):
        """
        Connect ASG and download designed sequences data into it
        :param asg_data: a List[List[float]] data type representing
        """
        is_connected = self._asg.connect()  # auto stop
        if is_connected == 1:
            self._asg.download_ASG_pulse_data(asg_data, [len(row) for row in asg_data])
        else:
            raise ConnectionError('ASG not connected')

    def download_asg_sequences(self, laser_seq: List[int] = None, mw_seq: List[int] = None,
                               tagger_seq: List[int] = None, N: int = 100000):
        """
        Download control sequences into the memory of ASG
        :param laser_seq: laser control sequence
        :param mw_seq: MW control sequence
        :param tagger_seq: tagger readout control sequence
        :param N: repetition number of sequences periods for each detection point
        """
        sequences = [laser_seq, mw_seq, tagger_seq]
        if not any(sequences):
            raise ValueError('laser_seq, mw_seq and tagger_seq cannot be all None')
        sequences = [seq for seq in sequences if seq is not None]  # non-None sequences
        lengths = np.unique(list(map(sum, sequences)))
        if len(lengths) != 1:
            raise ValueError('input sequences should have the same length')
        t = lengths.item()

        # configure ASG period information
        self._conf_time_paras(t, N)

        idx_laser_channel = self.channel['laser'] - 1
        idx_mw_channel = self.channel['mw'] - 1
        idx_tagger_channel = self.channel['tagger'] - 1

        self.reset_asg_sequence()
        if laser_seq is not None:
            self._asg_sequences[idx_laser_channel] = laser_seq
        if mw_seq is not None:
            self._asg_sequences[idx_mw_channel] = mw_seq
        if tagger_seq is not None:
            self._asg_sequences[idx_tagger_channel] = tagger_seq

        # connect & download pulse data
        self.asg_connect_and_download_data(self._asg_sequences)

    def configure_tagger_counting(self, apd_channel: int = None, asg_channel: int = None, reader: str = 'counter'):
        """
        Configure asg-channel and apd-channel for ASG. For Swabian Time Tagger, channel number range: [1, 8].
        :param apd_channel: APD channel number
        :param asg_channel: ASG channel number
        :param reader: counter of specific readout type
        """

        if apd_channel is not None:
            self.tagger_input['apd'] = apd_channel
        if asg_channel is not None:
            self.tagger_input['asg'] = asg_channel
        print('Current Tagger input channels:', self.tagger_input)

        # construct & execute Measurement instance
        N = self._asg_conf['N']

        if reader == 'counter':
            # continuous counting
            t_ps = int(self._asg_conf['t'] / C.pico)
            self.counter = tt.Counter(self.tagger, channels=[self.tagger_input['apd']], binwidth=t_ps, n_values=N)
        elif reader == 'cbm':
            # pulse readout
            if self.two_pulse_readout:
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

    def _start_device(self):
        """
        Start device: MW, ASG; Execute Measurement instance.
        """
        # 1. run ASG firstly
        self._data.clear()
        self._data_ref.clear()
        self._asg.start()

        # 2. restart self.counter if necessary
        if not self.counter.isRunning():
            self.counter.start()

        # 3. run MW then
        self._mw_instr.write_bool('OUTPUT:STATE', True)
        print('MW on/off status:', self._mw_instr.instrument_status_checking)

    def _get_data(self):
        self._data.append(self.counter.getData().ravel().tolist())
        self.counter.clear()

    def _get_data_ref(self):
        self._data_ref.append(self.counter.getData().ravel().tolist())
        self.counter.clear()

    def run(self):
        """
        Rough scheduling method
        1) start device
        2) acquire data timely
        """
        self._start_device()
        self._acquire_data()
        self.stop()

    def stop(self):
        """
        Stop hardware (ASG, MW, Tagger) scheduling
        """
        self.counter.stop()
        self._asg.stop()
        self._mw_instr.write_bool('OUTPUT:STATE', False)
        print('Stopped: Scheduling process has stopped')

    def close(self):
        """
        Release instrument (ASG, MW, Tagger) resources
        """
        self._asg.close_device()
        self._mw_instr.close()
        tt.freeTimeTagger(self.tagger)
        print('Closed: All instrument resources has been released')

    def configure_mw_paras(self, power: float = None, freq: float = None, regulate_pi: bool = False, *args, **kwargs):
        """
        Configure parameters of MW instrument
        :param power: power of MW, unit: dBm
        :param freq: frequency of MW, unit: Hz
        :param regulate_pi: whether regulate built-int MW pi pulse of this Scheduler
        """
        if power is not None:
            self._mw_conf['power'] = power
            self._mw_instr.write_float('POW', power)
            if regulate_pi:  # regulate time duration based on MW power
                self._regulate_pi_pulse(power=power)  # by power
        if freq is not None:
            self._mw_conf['freq'] = freq
            self._mw_instr.write_float('FREQUENCY', freq)
        self._mw_instr.write_bool('OUTPUT:STATE', True)

    def _regulate_pi_pulse(self, power: float = None, time: float = None):
        """
        Regulate time duration of MW pi pulse according to designed MW power, or vice verse
        :param power: MW power, unit: dBm
        :param time: time duration of MW pi pulse, unit: s
        :return:
        """
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
        """
        Reset all channels of ASG as ZERO signals
        """
        self._asg_sequences = [[0, 0] for _ in range(8)]

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

    def laser_on_seq(self):
        """
        Set sequence to control Laser keeping on during the whole period
        """
        idx_laser_channel = self.channel['laser'] - 1
        t = sum(self._asg_sequences[idx_laser_channel])
        self._asg_sequences[idx_laser_channel] = [t, 0]
        self.asg_connect_and_download_data(self._asg_sequences)
        self.asg.start()

    def laser_off_seq(self):
        """
        Set sequence to control Laser keeping off during the whole period
        """
        idx_laser_channel = self.channel['laser'] - 1
        self._asg_sequences[idx_laser_channel] = [0, 0]
        self.asg_connect_and_download_data(self._asg_sequences)
        self.asg.start()

    def mw_on_seq(self):
        """
        Set sequence to control MW keeping on during the whole period
        """
        idx_mw_channel = self.channel['mw'] - 1
        t = sum(self._asg_sequences[idx_mw_channel])
        if self.mw_ttl == 0:
            mw_seq = [0, t]
        else:
            mw_seq = [t, 0]
        self._asg_sequences[idx_mw_channel] = mw_seq
        self.asg_connect_and_download_data(self._asg_sequences)
        self.asg.start()

    def mw_off_seq(self):
        """
        Set sequence to control MW keeping off during the whole period
        """
        mw_seq = [0, 0]
        idx_mw_channel = self.channel['mw'] - 1
        self._asg_sequences[idx_mw_channel] = mw_seq
        self.asg_connect_and_download_data(self._asg_sequences)
        self.asg.start()

    def mw_control_seq(self, mw_seq: List[int] = None) -> Optional[List[int]]:
        """
        Get or set current MW control sequence
        :param mw_seq: designed MW control sequence, optional parameter
        :return: return current MW control sequence when mw_seq designed, otherwise return None
        """
        idx_mw_channel = self.channel['mw'] - 1
        if mw_seq is None:
            return self._asg_sequences[idx_mw_channel]
        else:
            self._asg_sequences[idx_mw_channel] = mw_seq
            self.asg_connect_and_download_data(self._asg_sequences)
            self._asg.start()

    def _conf_time_paras(self, t, N):
        """
        Configure characteristic time parameters
        :param t: total time of one ASG sequence period, unit: ns
        :param N: repetition number of ASG sequence
        """
        self._asg_conf['t'] = t * C.nano  # unit: s
        self._asg_conf['N'] = N
        self.asg_dwell = self._asg_conf['N'] * self._asg_conf['t']  # duration without padding
        self.time_pad = 0.035 * self.asg_dwell  # 0.035

    def _cal_counts_result(self):
        """
        Calculate counts of data, with respect to frequencies (scanning mode)
        If using two-pulse-readout strategy, calculate both signals and reference signals
        """
        if isinstance(self, FrequencyDomainScheduler):
            xs_name = 'freqs'
            xs = self._freqs
        elif isinstance(self, TimeDomainScheduler):
            xs_name = 'times'
            xs = self._times
        else:
            raise TypeError('unsupported function in this scheduler type')

        if self.with_ref:
            counts = [np.mean(ls) for ls in self._data]
            counts_ref = [np.mean(ls) for ls in self._data_ref]
            self._result = [xs, counts, counts_ref]
            self._result_detail = {
                xs_name: xs,
                'counts': counts,
                'counts_ref': counts_ref,
                'origin_data': self._data,
                'origin_data_ref': self._data_ref
            }
        else:
            if self.two_pulse_readout:
                counts_pairs = [(np.mean(ls[1::2]), np.mean(ls[::2])) for ls in self._data]
                counts = list(map(min, counts_pairs))
                counts_ref = list(map(max, counts_pairs))
                self._result = [xs, counts, counts_ref]
                self._result_detail = {
                    xs_name: xs,
                    'counts': counts,
                    'counts_ref': counts_ref,
                    'origin_data': self._data,
                }
            else:
                counts = [np.mean(ls) for ls in self._data]
                self._result = [xs, counts]
                self._result_detail = {
                    xs_name: xs,
                    'counts': counts,
                    'origin_data': self._data
                }

    def _gene_data_result_fname(self) -> str:
        """
        Generate file name of data acquisition result, based on time, data and random numbers
        :return: file name, str type
        """
        if self.with_ref:
            # calculate signals and reference counts
            fname = os.path.join(self.output_dir,
                                 '{}-counts-with-ref-{}-{}'.format(self.name.split()[0], str(datetime.date.today()),
                                                                   uuid.uuid1()))
        else:
            # just calculate signal counts
            fname = os.path.join(self.output_dir,
                                 '{}-counts-{}-{}'.format(self.name.split()[0], str(datetime.date.today()),
                                                          uuid.uuid1()))
        return fname

    def save_result(self, fname: str = None):
        """
        Save self.result_detail property
        :param fname: if not designed, will be randomly generated
        """
        if not self._result or not self._result_detail:
            raise ValueError('empty result cannot be saved')
        if fname is None:
            fname = os.path.join(self.output_dir,
                                 '{}-result-{}-{}'.format(self.name, str(datetime.date.today()), uuid.uuid1()))
        with open(fname + '.json', 'w') as f:
            json.dump(self._result_detail, f)
        print('Detailed data result has been saved into {}'.format(fname + '.json'))

    def __str__(self):
        return self.name

    @abc.abstractmethod
    def configure_odmr_seq(self, *args, **kwargs):
        """
        Configure ODMR detection sequences parameters
        """
        raise NotImplementedError

    @abc.abstractmethod
    def _acquire_data(self, *args, **kwargs):
        """
        Acquire data in real time
        """
        raise NotImplementedError

    @abc.abstractmethod
    def run_single_step(self, *args, **kwargs):
        """
        Run the Scheduler at one single detection parameter setting point
        """
        raise NotImplementedError

    @abc.abstractmethod
    def run_scanning(self, *args, **kwargs):
        """
        Scanning mode to run this Scheduler
        """
        raise NotImplementedError

    @property
    def laser_control(self) -> bool:
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
    def result(self) -> List[List[float]]:
        return self._result

    @property
    def result_detail(self) -> dict:
        return self._result_detail

    @property
    def sequences_strings(self) -> str:
        return seq_to_str(self._asg_sequences)

    @property
    def sequences_figure(self) -> Figure:
        return seq_to_fig(self._asg_sequences)


class FrequencyDomainScheduler(Scheduler):
    """
    Frequency-domain ODMR detection abstract class
    """

    def __init__(self, *args, **kwargs):
        super(FrequencyDomainScheduler, self).__init__(*args, **kwargs)
        self.name = 'Base ODMR Scheduler'

    def set_mw_freqs(self, start, end, step):
        """
        Set frequencies for scanning detection (e.g. CW-ODMR, Pulse-ODMR)
        All unit is "Hz"
        :param start: start frequency
        :param end: end frequency
        :param step: frequency step
        """
        # unit: Hz
        # n_freqs = int((end - start) / step + 1)
        self._freqs = np.arange(start, end + step / 2, step).tolist()
        n_freqs = len(self._freqs)
        if self.asg_dwell == 0:
            raise ValueError('"asg_dwell" is 0.0 currently. Please set ODMR sequences parameters firstly.')
        else:
            self.time_total = self.asg_dwell * n_freqs * 2 if self.with_ref else self.asg_dwell * n_freqs

    def _scan_freqs_and_get_data(self):
        """
        Scanning frequencies & getting data of Counter
        """
        # ===================================================================
        for _ in range(self.epoch_omit):
            self._mw_instr.write_float('FREQUENCY', self._freqs[0])
            # self._mw_instr.write_bool('OUTPUT:STATE', True)
            print('scanning freq {:.3f} GHz (trivial)'.format(self._freqs[0] / C.giga))
            time.sleep(self.time_pad + self.asg_dwell)

            if self.with_ref:
                # self._mw_instr.write_bool('OUTPUT:STATE', False)
                time.sleep(self.time_pad + self.asg_dwell)
        # ===================================================================

        mw_on_seq = self._asg_sequences[self.channel['mw'] - 1]
        for i, freq in enumerate(self._freqs):
            self._mw_instr.write_float('FREQUENCY', freq)
            # self._mw_instr.write_bool('OUTPUT:STATE', True)

            print('scanning freq {:.3f} GHz'.format(freq / C.giga))
            t = threading.Thread(target=self._get_data, name='thread-{}'.format(i))
            time.sleep(self.time_pad)
            time.sleep(self.asg_dwell)  # accumulate counts
            t.start()  # begin readout

            if self.with_ref:
                # modify the sequences
                self.mw_control_seq([0, 0])

                # reference data acquisition
                tr = threading.Thread(target=self._get_data_ref, name='thread-ref-{}'.format(i))
                time.sleep(self.time_pad)
                time.sleep(self.asg_dwell)
                tr.start()

                # recover the sequences
                self.mw_control_seq(mw_on_seq)

        print('finished data acquisition')

    def _acquire_data(self, *args, **kwargs):
        """
        Scanning time intervals to acquire data for Time-domain Scheduler
        :param with_ref: if True, MW 'on' and 'off' in turn
        """
        # 1. scan time intervals
        self._scan_freqs_and_get_data()

        # 2. calculate result (count with/without reference)
        self._cal_counts_result()

        # 3. save result
        self.save_result(self._gene_data_result_fname())

    def run_scanning(self, mw_control: str = 'on'):
        """
        Run the scheduler under scanning-frequency mode
        1) start device
        2) acquire data timely
        :param mw_control: 'on' or 'off'
        """
        mw_seq_on = self.mw_control_seq()
        if mw_control == 'off':
            self.mw_control_seq([0, 0])
        elif mw_control == 'on':
            pass
        else:
            raise ValueError('unsupported MW control parameter (should be "on" or "off"')

        print('Begin to run {}. Frequency: {:.3f} - {:.3f} GHz.'.format(self.name, self._freqs[0] / C.giga,
                                                                        self._freqs[-1] / C.giga))
        print('t: {:.2f} ns, N: {}, T: {:.2f} s, n_freqs: {}'.format(self._asg_conf['t'] / C.nano, self._asg_conf['N'],
                                                                     self.asg_dwell, len(self._freqs)))
        print('Estimated total running time: {:.2f} s'.format(self.time_total))

        self._start_device()
        self._acquire_data()  # scanning MW frequencies in this loop
        self.stop()

        # recover the asg control sequence for MW to be 'on'
        if mw_control == 'off':
            self.mw_control_seq(mw_seq_on)

    @abc.abstractmethod
    def configure_odmr_seq(self, *args, **kwargs):
        raise NotImplementedError

    def run_single_step(self, *args, **kwargs):
        raise NotImplementedError


class TimeDomainScheduler(Scheduler):
    """
    Time-domain ODMR detection abstract class
    """

    def __init__(self, *args, **kwargs):
        super(TimeDomainScheduler, self).__init__(*args, **kwargs)
        self.name = 'Time-domain ODMR Scheduler'

    def set_delay_times(self, start=None, end=None, step=None, times=None):
        """
        Set time intervals for scanning detection (e.g. Ramsey, Rabi, DD, Relaxation)
        All unit is "ns"
        :param start: start time interval
        :param end: end time interval
        :param step: time interval step
        :param times: time duration list
        """
        if times is not None:
            self._times = list(times)
        else:
            self._times = np.arange(start, end + step / 2, step).tolist()
        N = self._asg_conf['N']
        if N is None:
            raise ValueError('"N" is None currently. Please set ODMR sequences parameters firstly.')
        else:
            self.time_total = sum(self._times) * C.nano * N  # estimated total time

    def _scan_times_and_get_data(self):
        """
        Scanning time intervals & getting data of Counter
        """
        # =======================================================
        for _ in range(self.epoch_omit):
            # self._mw_instr.write_bool('OUTPUT:STATE', True)
            print('scanning freq {:.3f} ns (trivial)'.format(self._times[0]))
            self._gene_detect_seq(self._times[0])
            self._asg.start()
            time.sleep(self.time_pad + self.asg_dwell)

            if self.with_ref:
                # self._mw_instr.write_bool('OUTPUT:STATE', False)
                time.sleep(self.time_pad + self.asg_dwell)

        # =======================================================
        for i, duration in enumerate(self._times):
            self._gene_detect_seq(duration)
            self._asg.start()
            # self._mw_instr.write_bool('OUTPUT:STATE', True)
            print('scanning time interval: {:.3f} ns'.format(duration))

            # Signal readout
            t = threading.Thread(target=self._get_data, name='thread-{}'.format(i))
            time.sleep(self.time_pad)
            time.sleep(self.asg_dwell)  # accumulate counts
            t.start()  # begin readout

            # Reference readout
            if self.with_ref:
                self.mw_control_seq([0, 0])
                # self._mw_instr.write_bool('OUTPUT:STATE', False)
                tr = threading.Thread(target=self._get_data_ref, name='thread-ref-{}'.format(i))
                time.sleep(self.time_pad)
                time.sleep(self.asg_dwell)
                tr.start()

        print('finished data acquisition')

    def _acquire_data(self, *args, **kwargs):
        """
        Scanning time intervals to acquire data for Time-domain Scheduler
        """
        # 1. scan time intervals
        self._scan_times_and_get_data()

        # 2. calculate result (count with/without reference)
        print('two pulses:', self.two_pulse_readout)
        self._cal_counts_result()

        # 3. save result
        self.save_result(self._gene_data_result_fname())

    def run_scanning(self):
        """
        Run the scheduler under scanning-time-interval mode
        1) start device
        2) acquire data timely
        """
        print('Begin to run {}. Time intervals: {:.3f} - {:.3f} ns.'.format(self.name, self._times[0], self._times[-1]))
        print('N: {}, n_times: {}'.format(self._asg_conf['N'], len(self._times)))
        print('Estimated total running time: {:.2f} s'.format(self.time_total))

        self._start_device()
        self._acquire_data()  # scanning time intervals in this loop
        self.stop()

    def _gene_pseudo_detect_seq(self):
        """
        Generate pseudo pulses for visualization and regulation
        """
        ts = list(self._cache.values())[:-1]
        t_sum = sum([t for t in ts if t is not None])
        self._gene_detect_seq(int(t_sum / 40) * 10)

    @abc.abstractmethod
    def _gene_detect_seq(self, *args, **kwargs):
        raise NotImplementedError

    @abc.abstractmethod
    def configure_odmr_seq(self, *args, **kwargs):
        raise NotImplementedError

    def run_single_step(self, t_free) -> List[float]:
        """
        Single-interval & single-power setting for running the scheduler once
        :param t_free: free precession time, unit: ns
        :return: 1-D array-like data: [N,]
        """
        print('running with time = {:.3f} ns, MW power = {:.2f} dBm ...'.format(self.asg_dwell, self._mw_conf['power']))

        # generate ASG sequences
        self._gene_detect_seq(t_free)

        # start sequence for time: N * t
        self._start_device()
        time.sleep(2)  # 先让激光和微波开几秒
        time.sleep(self.asg_dwell)
        counts = self.counter.getData().ravel().tolist()
        self.stop()

        return counts
