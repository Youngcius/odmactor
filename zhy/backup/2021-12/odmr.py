"""
Backup on 2021-11-23
"""
import abc
from odmactor.scheduler.base import Scheduler
import datetime
import time
import os
import numpy as np
import scipy.constants as C
import TimeTagger as tt
import threading
import pickle


class ODMRScheduler(Scheduler):
    """
    ODMR measurement abstract class
    """

    def __init__(self, *args, **kwargs):
        super(ODMRScheduler, self).__init__(*args, **kwargs)
        self.name = 'Base ODMR Scheduler'

    def _start_device(self):
        # 1. run ASG firstly
        self._data.clear()
        self._asg.start()

        # restart self.counter if necessary
        if not self.counter.isRunning():
            self.counter.start()

        # 2. run MW then
        self._mw_instr.write_bool('OUTPUT:STATE', True)
        # if self.mw_exec_mode == 'scan-center-span' or self.mw_exec_mode == 'scan-start-stop':
        #     self._mw_instr.write_str('SWE:FREQ:EXEC')  # trigger the sweep
        #     raise warnings.warn('Using the SWEEP mode of MW', DeprecationWarning)
        #       self._mw_instr.write_str('SOUR:PULM:TRIG:MODE SING')
        #       self._mw_instr.write_str('SOUR:PULM:TRIG:IMM')
        print('MW on/off status:', self._mw_instr.instrument_status_checking)

    def _scan_freqs_and_get_data(self):
        """
        Sanning frequencies & getting data of Counter
        """
        for i, freq in enumerate(self._freqs):
            self._mw_instr.write_float('FREQUENCY', freq)
            print('scanning freq {:.4f} GHz'.format(freq / C.giga))
            t = threading.Thread(target=self._get_data, name='thread-{}'.format(i))
            time.sleep(self.time_pad)
            time.sleep(self.asg_dwell)  # accumulate counts
            t.start()  # begin readout
            time.sleep(self.time_pad)
            t.join()
        print('finished data acquisition')

    def run_scanning(self, mw_control='on'):
        """
        Run the scheduler under scanning-frequency mode
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


class CWScheduler(ODMRScheduler):
    """
    Continuous-wave ODMR scheduler
    """

    def __init__(self, *args, **kwargs):
        super(CWScheduler, self).__init__(*args, **kwargs)
        self.name = 'CW ODMR Scheduler'

    # def configure_odmr_seq(self, t, t_read, N: int):
    #     """
    #     Wave form for single period:
    #         laser (no asg control sequence):
    #         ---------------------------------
    #         |                               |
    #         |                               |
    #         microwave (no asg control sequence):
    #         ---------------------------------
    #         |                               |
    #         |                               |
    #         asg tagger acquisition channel:
    #                                       ---
    #                                       | |
    #         ------------------------------| |
    #     All units for the parameters is 'ns'
    #     :param t: time for a single period of ASG sequence
    #     :param t_read: time span for fluorescence signal readout
    #     :param N: number of ASG operation periods
    #     """
    #     # unit: ns
    #     self._asg_conf['t'] = t * C.nano  # ns --> s
    #     self._asg_conf['N'] = N
    #     self.asg_dwell = self._asg_conf['N'] * self._asg_conf['t']
    #     self.mw_dwell = self.asg_dwell + self.time_pad*2
    #     # generate ASG wave forms
    #     idx_tagger_channel = self.channel['tagger'] - 1
    #     idx_apd_channel = self.channel['apd'] - 1
    #     tagger_seq = [0, t - t_read, t_read, 0]
    #     self._asg_sequences[idx_tagger_channel] = tagger_seq
    #     self._asg_sequences[idx_apd_channel] = tagger_seq
    #
    #     # connect & download pulse data
    #     self.asg_connect_and_download_data(self._asg_sequences)

    def configure_odmr_seq(self, t, N: int):
        """
        Wave form for single period:
            laser (no asg control sequence):
            ---------------------------------
            |                               |
            |                               |
            microwave (no asg control sequence):
            ---------------------------------
            |                               |
            |                               |
            asg tagger acquisition channel (continuous readout in fact):
            ---------------------------------
            |                               |
            |                               |
        All units for the parameters is 'ns'
        :param t: binwidth parameter for TimeTagger.Counter
        :param N: n_values for TimeTagger.Counter
        """
        # unit: ns
        self._conf_time_paras(t, N)

        # generate ASG wave forms ('pulse' mode for Laser)
        idx_laser_channel = self.channel['laser'] - 1
        idx_mw_channel = self.channel['mw'] - 1
        idx_apd_channel = self.channel['apd'] - 1
        continuous_seq = [t, 0]
        self._asg_sequences = [[0, 0] for i in range(8)]

        self._asg_sequences[idx_laser_channel] = continuous_seq
        self._asg_sequences[idx_mw_channel] = continuous_seq
        self._asg_sequences[idx_apd_channel] = continuous_seq

        # connect & download pulse data
        self.asg_connect_and_download_data(self._asg_sequences)

    def _acquire_data(self, *args, **kwargs):
        # 1. scan freq
        print('CW _acquire_data')
        self._scan_freqs_and_get_data()
        # 2. calculate result
        self._cal_counts_result()
        # 3. save result
        self.save_result()

    # def _acquire_data(self, *args, **kwargs):
    #
    #     # for i, freq in enumerate(self._freqs):
    #     #     print('scanning freq {:.4f} GHz'.format(freq / C.giga))
    #     n_freqs = len(self._freqs)
    #     beg = time.time_ns()
    #     for i in range(n_freqs):
    #         t = threading.Thread(target=self._get_data, name='thread-{}'.format(i))
    #         time.sleep(self.time_pad)
    #         time.sleep(self.asg_dwell)  # accumulate counts
    #         t.start()  # begin readout
    #         time.sleep(self.time_pad)
    #         t.join()
    #     end = time.time_ns()
    #     self.sync_delay = end - beg
    #     print('ideal time', self.mw_dwell * len(self._freqs))
    #     print('actualt time', self.sync_delay)
    #
    #     # 计算计数
    #     counts = [np.mean(ls) for ls in self._data]
    #     self._result = [self._freqs, counts]
    #     self._result_detail = {
    #         'freqs': self._freqs,
    #         'counts': counts,
    #         'origtin_data': self._data
    #     }
    #
    #     fname = os.path.join(self.output_dir,
    #                          'CW-ODMR-result-{}-{}'.format(str(datetime.date.today()), round(time.time() / 120)))
    #     with open(fname + '.json', 'w') as f:
    #         json.dump(self._result_detail, f)
    #     np.savetxt(fname + '.txt', self._result)
    #     print('result has been saved into {}'.format(fname + '.json'))

    # def run(self, mw_control='on'):
    #     """
    #     1) start device
    #     2) acquire data timely
    #     """
    #     mw_seq_on = self._asg_sequences[self.channel['mw'] - 1]
    #     if mw_control == 'off':
    #         # mw_seq_off = [0, sum(mw_seq_on)]
    #         self._asg_sequences[self.channel['mw'] - 1] = [0, 0]
    #         self.asg_connect_and_download_data(self._asg_sequences)
    #     elif mw_control == 'on':
    #         pass
    #     else:
    #         raise ValueError('unsupported MW control parameter (should be "or" or "off"')
    #     print('Begin to run {}. Frequency: {:.4f} - {:.4f} GHz.'.format(self.name, self._freqs[0], self._freqs[-1]))
    #     self._start_device()
    #     self._acquire_data()
    #     self.stop()
    #
    #     # 恢复微波的ASG的MW通道为 on
    #     self._asg_sequences[self.channel['mw'] - 1] = mw_seq_on
    #     self.asg_connect_and_download_data(self._asg_sequences)

    def run_origin(self, mw_control='on'):
        """
        1) start device
        2) acquire data timely
        """
        mw_seq_on = self._asg_sequences[self.channel['mw'] - 1]
        if mw_control == 'off':
            mw_seq_off = [0, sum(mw_seq_on)]
            self._asg_sequences[self.channel['mw'] - 1] = mw_seq_off
            self.asg_connect_and_download_data(self._asg_sequences)
        elif mw_control == 'on':
            pass
        else:
            raise ValueError('unsupported MW control parameter (should be "or" or "off"')

        self._start_device()
        self._acquire_data()

        # 恢复微波的ASG的MW通道为 on
        self._asg_sequences[self.channel['mw'] - 1] = mw_seq_on
        self.asg_connect_and_download_data(self._asg_sequences)

        # self._start_device()
        # self._acquire_data()

        data_sig = []
        data_ref = []
        t = self._asg_conf['t']  # unit: s
        num_freqs = len(self._freqs)
        self.counter = tt.Counter(self.tagger, channels=[1], binwidth=int(self._asg_conf['t'] / C.pico),
                                  n_values=self._asg_conf['N'])

        for freq in self._freqs:
            print('scanning freq {:.4f} GHz'.format(freq / C.giga))
            # self._mw_instr.write_bool('OUTPUT:STATE', True)
            self._mw_instr.write_float('FREQUENCY', freq)
            time.sleep(0.2)

            # self.laser_on_seq(t / C.nano)
            self.mw_on_seq(t / C.nano)
            time.sleep(self.asg_dwell)
            data_sig.append(self.counter.getData())
            self.counter.clear()
            # print('---------- on')
            # print(self.counter.getData().ravel())
            # print(self._asg_sequences)
            # self.counter.clear()
            # print('after clear:', self.counter.getData())

            self.asg.stop()  # stop firstly
            # self.laser_off_seq()
            self.mw_off_seq()
            time.sleep(self.asg_dwell)
            data_ref.append(self.counter.getData())
            self.counter.clear()
            # print('---------- off')
            # print(self.counter.getData().ravel())
            # print(self._asg_sequences)
            # self.counter.clear()
            # print('after clear:', self.counter.getData())

        self.asg.stop()
        self.counter.stop()

        data_sig_avg = [np.mean(ls) for ls in data_sig]
        data_ref_avg = [np.mean(ls) for ls in data_ref]
        contrast = [np.abs(data_sig_avg[i] - data_ref_avg[i]) / data_ref_avg[i] for i in range(num_freqs)]

        self._result = np.array([self._freqs, contrast])

        self._result_detail = {
            'freqs': self._freqs,
            'contrast': contrast,
            'data_ref': data_sig,
            'data_sig': data_ref
        }

        fname = os.path.join(self.output_dir,
                             'CW-ODMR-result-{}-{}'.format(str(datetime.date.today()), round(time.time())))
        with open(fname + '.pkl', 'wb') as f:
            pickle.dump(self._result_detail, f)
        np.savetxt(fname + '.txt', self._result)
        print('data has been saved into {}'.format(fname))

    def run_single_step(self, power, freq, mw_control='on'):
        """
        Single-frequency & single-power setting for running the scheduler
        :param power: MW power, unit: dBm
        :param freq: MW frequency, unit: Hz
        :param mw_control: 'on' or 'off'
        :return: 1-D array: [N,]
        """
        print('running for freq = {:.4f} GHz, power = {:.2f} dBm ...'.format(freq / C.giga, power))

        # set MW parameters
        self.configure_mw_paras(power, freq)

        # start sequence for time: N*t
        if mw_control == 'on':
            self.mw_on_seq(self._asg_conf['t'] / C.nano)  # auto-start ASG
        elif mw_control == 'off':
            self.mw_off_seq()  # auto-start ASG
        else:
            raise ValueError('unsupported mw_control parameter')

        time.sleep(3)  # 先让激光和微波开几秒
        self.counter.start()
        time.sleep(self.time_pad)
        time.sleep(self.asg_dwell)
        data = self.counter.getData().ravel()
        self.stop()

        return data

    def mw_on_seq(self, t):

        mw_seq = [t, 0]
        idx_mw_channel = self.channel['mw'] - 1
        self._asg_sequences[idx_mw_channel] = mw_seq
        # self.asg.close_device()

        self.asg_connect_and_download_data(self._asg_sequences)

        self.asg.start()

    def mw_off_seq(self):

        mw_seq = [0, 0]
        idx_mw_channel = self.channel['mw'] - 1
        self._asg_sequences[idx_mw_channel] = mw_seq
        # self.asg.close_device()
        self.asg_connect_and_download_data(self._asg_sequences)
        self.asg.start()


class PulseScheduler(ODMRScheduler):
    """
    Pulse-based ODMR manipulation scheduler
    """

    def __init__(self, *args, **kwargs):
        super(PulseScheduler, self).__init__(*args, **kwargs)
        self.name = 'Pulse ODMR Scheduler'
        self.two_pulse_readout = False
        # self.transition_time = transition_time  # 计数过渡到平衡时的时间（开关微波时）

    def configure_odmr_seq(self, t_init, t_mw, t_read_sig, t_read_ref=None, inter_init_mw=3000,
                           inter_readout=200, inter_period=200, N: int = 1000):
        """
        Wave form for single period:
            asg laser channel:
            -----               ---------
            |   |               |       |
            |   |---------------|       |----
            asg microwave channel:
                    -------------
                    |           |
            --------|           |------------
            asg tagger acquisition channel:
                                ---   ---
                                | |   | |
            --------------------| |---| |----
        All units for the parameters is 'ns'
        :param t_init: time span for laser initialization, e.g. 5000
        :param t_mw: time span for microwave actual operation in a ASG period, e.g. 800
        :param t_read_sig: time span for fluorescence signal readout, e.g. 400
        :param t_read_ref: time span for reference signal readout, e.g. 400
                            if the parameter is not assigned, means single-pulse readout
        :param inter_init_mw: time interval between laser initialization and MW operation pulses, e.g. 3000
        :param inter_readout: interval between single readout pulse and reference signal readout, e.g. 200
                            when t_read_ref is assigned, this parameter will play its role
        :param inter_period: interval between two neighbor periods, e.g. 200
        :param N: number of ASG operation periods
        """
        # unit: ns
        # total time for 'N' period, also for MW operation time at each frequency point
        if t_read_ref is not None:
            self.two_pulse_readout = True
            t = t_init + inter_init_mw + t_mw + t_read_sig + inter_readout + t_read_ref + inter_period
            laser_seq = [t_init, inter_init_mw + t_mw, t_read_sig + inter_readout + t_read_ref, inter_period]
            mw_seq = [0, t_init + inter_init_mw, t_mw, t_read_sig + inter_readout + t_read_ref + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw, t_read_sig, inter_readout, t_read_ref, inter_period]
        else:
            # single-pulse readout
            t = t_init + inter_init_mw + t_mw + t_read_sig + inter_period
            laser_seq = [t_init, inter_init_mw + t_mw, t_read_sig, inter_period]
            mw_seq = [0, t_init + inter_init_mw, t_mw, t_read_sig + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw, t_read_sig, inter_period]

        self._conf_time_paras(t, N)

        # generate ASG wave forms
        idx_laser_channel = self.channel['laser'] - 1
        idx_mw_channel = self.channel['mw'] - 1
        idx_tagger_channel = self.channel['tagger'] - 1
        idx_apd_channel = self.channel['apd'] - 1

        self._asg_sequences = [[0, 0] for i in range(8)]
        self._asg_sequences[idx_laser_channel] = laser_seq
        self._asg_sequences[idx_mw_channel] = mw_seq
        self._asg_sequences[idx_tagger_channel] = tagger_seq
        self._asg_sequences[idx_apd_channel] = [t, 0]

        # connect & download pulse data
        self.asg_connect_and_download_data(self._asg_sequences)

    def _acquire_data(self):
        """
        Default setting: save file into text files.
        """
        # 1. scan freq
        self._scan_freqs_and_get_data()
        # 2. calculate result (count or contrast)
        if self.two_pulse_readout:
            # calculate contrasts
            self._cal_contrasts_result()
            fname = os.path.join(self.output_dir,
                                 'Pulse-ODMR-contrasts-{}-{}'.format(str(datetime.date.today()),
                                                                     round(time.time())))
        else:
            # just calculate counts
            self._cal_counts_result()
            fname = os.path.join(self.output_dir,
                                 'Pulse-ODMR-counts-{}-{}'.format(str(datetime.date.today()), round(time.time())))
        # 3. save result
        self.save_result(fname)

    # def run(self, mw_control='on'):
    #     mw_seq_on = self._asg_sequences[self.channel['mw'] - 1]
    #     if mw_control == 'off':
    #         mw_seq_off = [0, sum(mw_seq_on)]
    #         self._asg_sequences[self.channel['mw'] - 1] = mw_seq_off
    #         self.asg_connect_and_download_data(self._asg_sequences)
    #     elif mw_control == 'on':
    #         pass
    #     else:
    #         pass
    #     print('Begin to run {}. Frequency: {:.4f} - {:.4f} GHz.'.format(self.name, self._freqs[0], self._freqs[-1]))
    #     print('Estimated total running time: {:.2f} s'.format(self.time_total))
    #     self._start_device()
    #     self._acquire_data()
    #     self.stop()
    #     # 恢复微波的ASG的MW通道为 on
    #     self._asg_sequences[self.channel['mw'] - 1] = mw_seq_on
    #     self.asg_connect_and_download_data(self._asg_sequences)
