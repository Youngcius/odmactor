from odmactor.scheduler.base import Scheduler
import TimeTagger as tt
import numpy as np
import os
import datetime
import time
import scipy.constants as C
from odmactor.utils import cut_edge_zeros, cal_contrast
import pickle

"""
Pulse detection (frequency-domain method)
输入：
	- 微波频率范围和步长（起止点/中心和宽度）--> freqs = linspace(,,)
	- ASG 序列周期个数，serial time intervals parameters. e.g. 10ms * 100 = 1s ，后者作为微波在单个频率点的持续时间 T = N * t
	- 单个周期内激光初始化时间 t_init, 微波翻转时间 t_mw，读出时间 t_read
	-
输出：对比度数据，
调度过程  结果读出 ：
	- 生成两通道数据（channel['laser']，channel['mw'] --> asg_data: Matrix）
	- For each freq in freqs: Asg start --> Asg stop 历时 T ，存储 APD 数据 到 self._cache
	- 最后统一计算各频率点的对比度，结果到 self.result : {'freqs': …, 'ratio': …}，self.result_detail
"""


class PulseScheduler(Scheduler):
    """
    Pulse-based ODMR manipulation scheduler
    """

    def __init__(self, *args, **kwargs):
        super(PulseScheduler, self).__init__(*args, **kwargs)
        self.name = 'Pulse ODMR Scheduler'

    def configure_odmr_seq(self, t_mw, t_init, t_read_sig, t_read_ref, t_interval=20, N: int = 100):
        """
        Wave form for single period:
            asg laser channel:
            -----                   ---------
            |   |                   |       |
            |   |-------------------|       |----
            asg microwave channel:
                    -------------
                    |           |
            --------|           |----------------
            asg tagger acquisition channel:
                                    ---   ---
                                    | |   | |
            ------------------------| |---| |----
        All units for the parameters is 'ns'
        :param t_mw: time span for microwave actual operation in a ASG period
        :param t_init: time span for laser initialization
        :param t_read_sig: time span for fluorescence signal readout
        :param t_read_ref: time span for reference signal readout
        :param t_interval: time span for interval
        :param N: number of ASG operation periods
        """
        # unit: ns
        t = t_init + t_interval + t_mw + t_interval + t_read_sig + t_interval + t_read_ref + t_interval
        # total time for 'N' period, also for MW operation time at each frequency point
        self._asg_conf['t'] = t * C.nano  # unit: s
        self._asg_conf['N'] = N
        self.asg_dwell = self._asg_conf['N'] * self._asg_conf['t']  # duration without padding
        self.mw_dwell = self.asg_dwell + self.time_pad

        # generate ASG wave forms
        idx_laser_channel = self.channel['laser'] - 1
        idx_mw_channel = self.channel['mw'] - 1
        idx_tagger_channel = self.channel['tagger'] - 1
        idx_apd_channel = self.channel['apd'] - 1
        laser_seq = [t_init, t_mw + 2 * t_interval, t_read_sig + t_interval + t_read_ref, t_interval]
        if self.mw_ttl == 1:
            mw_seq = [0, t_init + t_interval, t_mw, 3 * t_interval + t_read_sig + t_read_ref]
        else:
            mw_seq = [t_init + t_interval, t_mw, 3 * t_interval + t_read_sig + t_read_ref, 0]
        tagger_seq = [0, t_mw + 2 * t_interval + t_init, t_read_sig, t_interval, t_read_ref, t_interval]
        self._asg_sequences = [[0, 0] for i in range(8)]
        self._asg_sequences[idx_laser_channel] = laser_seq
        self._asg_sequences[idx_mw_channel] = mw_seq
        self._asg_sequences[idx_tagger_channel] = tagger_seq
        self._asg_sequences[idx_apd_channel] = tagger_seq

        # connect & download pulse data
        self.asg_connect_and_download_data(self._asg_sequences)

    def _start_device(self):
        # 1. run MW firstly
        self._mw_instr.write_bool('OUTPUT:STATE', True)
        if self.mw_exec_mode == 'scan-center-span' or self.mw_exec_mode == 'scan-start-stop':
            self._mw_instr.write_str('SWE:FREQ:EXEC')  # trigger the sweep
            # self._mw_instr.write_str('SOUR:PULM:TRIG:MODE SING')
            # self._mw_instr.write_str('SOUR:PULM:TRIG:IMM')

        beg = time.time_ns()

        # 2. run ASG then
        self._asg.start()

        end = time.time_ns()
        self.sync_delay = end - beg

        print('MW status now:', self._mw_instr.instrument_status_checking)

        # execute Measurement instance
        self.counter = tt.CountBetweenMarkers(self.tagger, self.tagger_input['apd'],
                                              begin_channel=self.tagger_input['asg'],
                                              end_channel=-self.tagger_input['asg'],
                                              n_values=self._asg_conf['N'] * 2)
        # self.counter.startFor(int((self.time_total + 5) / C.pico))  # parameter unit: ps

    def _acquire_data(self):
        """
        Default setting: save file into text files. TODO
        """

        # acquire data from each period, from TimeTagger TODO
        # 0~T: N*2 data points, N*2 < n_values
        # or: automatically store
        # T = self._asg_conf['t'] * self._asg_conf['N']
        # t_cur = 0
        # while t_cur < self.time_total:
        #     t_cur += T
        #     counts = self.counter.getData()
        #

        data = []
        # while self.counter.isRunning():
        #     time.sleep(self.mw_dwell)  # for each MW frequency
        #     data.append(self.counter.getData())
        #     self.counter.clear()

        for freq in self._freqs:
            print('scanning freq {:.2f} GHz'.format(freq / C.giga))
            self.counter.clear()
            time.sleep(self.time_pad / 2)
            time.sleep(self.asg_dwell)
            data.append(self.counter.getData())

            time.sleep(self.time_pad / 2)

        self.counter.stop()
        # contrast_list = [cal_contrast(ls) for ls in data]
        #
        # contrast = [np.mean(item) for item in contrast_list]  # average contrast
        contrast = [cal_contrast(ls) for ls in data]
        print(len(self._freqs), len(contrast))
        self._result = [self._freqs, contrast]
        self._result_detail = {
            'freqs': self._freqs,
            'contrast': contrast,
            'origin_data': data,
            # 'contrast_list': contrast_list
        }

        fname = os.path.join(self.output_dir,
                             'Pulse-ODMR-result-{}-{}'.format(str(datetime.date.today()), round(time.time() / 120)))
        with open(fname + '.pkl', 'wb') as f:
            pickle.dump(self._result_detail, f)
        # np.savetxt(fname + '.txt', self._result)
        print('data has been saved into {}'.format(fname))
