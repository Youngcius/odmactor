from odmactor.scheduler.base import Scheduler
import TimeTagger as tt
import numpy as np
import threading
import os
import json
import datetime
import time
import scipy.constants as C
from odmactor.utils import cut_edge_zeros, cal_contrast
import pickle

# 20u
# rise: 20s
# down: 20s
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
经验参数：
    - init: 5us
    - init-mw-interval 3us
    - mw: 5us
    - readout sign: 400ns
    - readout-interval: 200ns
    - readout ref: 400s
"""


class PulseScheduler(Scheduler):
    """
    Pulse-based ODMR manipulation scheduler
    """

    def __init__(self, transition_time=0, *args, **kwargs):
        super(PulseScheduler, self).__init__(*args, **kwargs)
        self.name = 'Pulse ODMR Scheduler'
        self.two_pulse_readout = False
        self.transition_time = transition_time  # 计数过渡到平衡时的时间（开关微波时）

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
        :param t_mw: time span for microwave actual operation in a ASG period, e.g. 5000
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

        self._asg_conf['t'] = t * C.nano  # unit: s
        self._asg_conf['N'] = N
        self.asg_dwell = self._asg_conf['N'] * self._asg_conf['t']  # duration without padding
        self.mw_dwell = self.asg_dwell + self.time_pad * 2 + self.transition_time

        # generate ASG wave forms
        idx_laser_channel = self.channel['laser'] - 1
        idx_mw_channel = self.channel['mw'] - 1
        idx_tagger_channel = self.channel['tagger'] - 1
        idx_apd_channel = self.channel['apd'] - 1

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
        print('MW status now:', self._mw_instr.instrument_status_checking)

        # beg = time.time_ns()

        # 2. run ASG then
        self._data.clear()
        self._asg.start()

        # end = time.time_ns()
        # self.sync_delay = end - beg

    def _acquire_data(self):
        """
        Default setting: save file into text files.
        """
        for i, freq in enumerate(self._freqs):
            # print('scanning freq {:.2f} GHz'.format(freq / C.giga))
            t = threading.Thread(target=self._get_data, name='readout-{}'.format(i))
            time.sleep(self.time_pad + self.transition_time)
            time.sleep(self.asg_dwell)  # effective time for accumulate counts
            t.start()  # another thread for readout
            time.sleep(self.time_pad)
            t.join()

        if self.two_pulse_readout:
            # calculate the contrasts
            self._cal_contrasts_result()
            fname = os.path.join(self.output_dir,
                                 'Pulse-ODMR-contrasts-{}-{}'.format(str(datetime.date.today()),
                                                                     round(time.time() / 120)))
        else:
            # just calculate the counts
            self._cal_counts_result()
            fname = os.path.join(self.output_dir,
                                 'Pulse-ODMR-counts-{}-{}'.format(str(datetime.date.today()), round(time.time() / 120)))
        # from pprint import pprint
        # pprint(self.result_detail)
        print(self.result)
        with open(fname + '.json', 'w') as f:
            json.dump(self._result_detail, f)
        np.savetxt(fname + '.txt', self._result)
        print('result has been saved into {}'.format(fname + '.json'))

    def run(self, mw_control='on'):
        mw_seq_on = self._asg_sequences[self.channel['mw'] - 1]
        if mw_control == 'off':
            mw_seq_off = [0, sum(mw_seq_on)]
            self._asg_sequences[self.channel['mw'] - 1] = mw_seq_off
            self.asg_connect_and_download_data(self._asg_sequences)
        elif mw_control == 'on':
            pass
        else:
            pass
        self._start_device()
        self._acquire_data()
        self.stop()
        # 恢复微波的ASG的MW通道为 on
        self._asg_sequences[self.channel['mw'] - 1] = mw_seq_on
        self._asg.stop()
        self.asg_connect_and_download_data(self._asg_sequences)

    def _cal_contrasts_result(self):
        contrasts = [cal_contrast(ls) for ls in self._data]
        self._result = [self._freqs, contrasts]
        self._result_detail = {
            'freqs': self._freqs,
            'contrast': contrasts,
            'origin_data': self._data,
        }

    def _cal_counts_result(self):
        counts = [np.mean(ls) for ls in self._data]
        self._result = [self._freqs, counts]
        self._result_detail = {
            'freqs': self._freqs,
            'counts': counts,
            'origin_data': self._data
        }
