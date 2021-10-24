import datetime
import time
import os
from odmactor.scheduler import Scheduler
import numpy as np
import scipy.constants as C
import TimeTagger as tt
import pickle

"""
Continuous-wave detection (frequency-domain method)
Basic: 激光和微波同时施加，单光子探测器持续探测全过 程中的光子数，通过扫描微波频 率产生频谱。
Remarks: 连续波谱不是典型 的量子传感过程，实验中系统处 于开放稳态，相干性在此实验中几乎不起作用
"""


class CWScheduler(Scheduler):
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
    #     self.mw_dwell = self.asg_dwell + self.time_pad
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
        self._asg_conf['t'] = t * C.nano  # ns --> s
        self._asg_conf['N'] = N
        self.asg_dwell = self._asg_conf['N'] * self._asg_conf['t']
        # self.mw_dwell = self.asg_dwell + self.time_pad
        self.mw_dwell = self.asg_dwell

        # generate ASG wave forms ('pulse' mode for Laser)
        idx_laser_channel = self.channel['laser'] - 1
        idx_mw_channel = self.channel['mw'] - 1
        # idx_tagger_channel = self.channel['tagger'] - 1
        idx_apd_channel = self.channel['apd'] - 1
        continuous_seq = [t, 0]
        self._asg_sequences = [[0, 0] for i in range(8)]

        self._asg_sequences[idx_laser_channel] = continuous_seq
        self._asg_sequences[idx_mw_channel] = continuous_seq
        # self._asg_sequences[idx_tagger_channel] = continuous_seq
        self._asg_sequences[idx_apd_channel] = continuous_seq

        # connect & download pulse data
        # self.asg_connect_and_download_data(self._asg_sequences)

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

    def _start_device(self, *args, **kwargs):

        # execute Measurement instance
        # self.counter = tt.CountBetweenMarkers(self.tagger, self.tagger_input['apd'],
        #                                       begin_channel=self.tagger_input['asg'],
        #                                       end_channel=-self.tagger_input['asg'],
        #                                       n_values=self._asg_conf['N'])
        # self.counter.startFor(int((self.time_total + 5) / C.pico))  # parameter unit: ps
        self.counter = tt.Counter(self.tagger, channels=[1], binwidth=int(self._asg_conf['t'] / C.pico),
                                  n_values=self._asg_conf['N'])

        # run MW firstly
        self._mw_instr.write_bool('OUTPUT:STATE', True)
        if self.mw_exec_mode == 'scan-center-span' or self.mw_exec_mode == 'scan-start-stop':
            self._mw_instr.write_str('SWE:FREQ:EXEC')  # trigger the sweep

        a = time.time_ns()
        # run ASG then
        self._asg.start()
        b = time.time_ns()
        self.sync_delay = b - a

    def _acquire_data(self, *args, **kwargs):
        data = []
        self.time_log = []
        # while self.counter.isRunning():
        for freq in self._freqs:
            # print('scanning freq {:.2f} GHz'.format(freq / C.giga))
            time.sleep(self.time_pad / 2)
            time.sleep(self.mw_dwell)  # for each MW frequency
            data.append(self.counter.getData())
            self.time_log.append(time.time_ns())
            # self.counter.clear()
            time.sleep(self.time_pad / 2)
        self.counter.stop()

        data_avg = [np.mean(ls) for ls in data]
        contrast = [item / data_avg[0] for item in data_avg]  # contrast reference: the first frequency point
        self._result = np.array([self._freqs, contrast])
        # self._result = [self._freqs, contrast]

        self._result_detail = {
            'freqs': self._freqs,
            'contrast': contrast,
            'origin_data': data
        }

        fname = os.path.join(self.output_dir,
                             'CW-ODMR-result-{}-{}'.format(str(datetime.date.today()), round(time.time() / 120)))
        with open(fname + '.pkl', 'wb') as f:
            pickle.dump(self._result_detail, f)
        np.savetxt(fname + '.txt', self._result)
        print('data has been saved into {}'.format(fname))

        np.savetxt('../output/time_log.txt', self.time_log)

    def run(self):
        """
        1) start device
        2) acquire data timely
        """
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
                             'CW-ODMR-result-{}-{}'.format(str(datetime.date.today()), round(time.time() / 120)))
        with open(fname + '.pkl', 'wb') as f:
            pickle.dump(self._result_detail, f)
        np.savetxt(fname + '.txt', self._result)
        print('data has been saved into {}'.format(fname))
