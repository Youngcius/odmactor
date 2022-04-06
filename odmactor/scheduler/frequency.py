"""
Frequency-domain ODMR detecting classes
---
1. Continuous-wave detection (frequency-domain method)
Basic: 激光和微波同时施加，单光子探测器持续探测全过程中的光子数，通过扫描微波频率产生频谱
Remarks: 连续波谱不是典型的量子传感过程，实验中系统处于开放稳态，相干性在此实验中几乎不起作用
---
2. Pulse detection (frequency-domain method)
"""

import threading
import time
import numpy as np
import scipy.constants as C
from typing import List
from odmactor.scheduler.base import FrequencyDomainScheduler
from odmactor import utils


class CWScheduler(FrequencyDomainScheduler):
    """
    Continuous-Wave ODMR scheduler
    """

    def __init__(self, *args, **kwargs):
        super(CWScheduler, self).__init__(*args, **kwargs)
        self.name = 'CW ODMR Scheduler'

    def configure_odmr_seq(self, period, N: int):
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
        if self.use_lockin:
            mw_seq = [period, period]
            cont_seq = [period * 2, 0]
            self.download_asg_sequences(laser_seq=cont_seq, mw_seq=mw_seq, lockin_seq=mw_seq, N=N)

        else:
            cont_seq = [period, 0]
            if self.mw_ttl == 0:
                mw_seq = utils.flip_sequence(cont_seq)
            else:
                mw_seq = cont_seq
            self.download_asg_sequences(laser_seq=cont_seq, mw_seq=mw_seq, tagger_seq=cont_seq, N=N)

    def run_single_step(self, power, freq, mw_control='on') -> List[float]:
        """
        Single-frequency & single-power setting for running the scheduler
        :param power: MW power, unit: dBm
        :param freq: MW frequency, unit: Hz
        :param mw_control: 'on' or 'off'
        :return: 1-D array: [N,]
        """
        print('running for freq = {:.3f} GHz, power = {:.2f} dBm ...'.format(freq / C.giga, power))

        # set MW parameters
        self.configure_mw_paras(power, freq)

        mw_seq_on = self.mw_control_seq()
        if mw_control == 'off':
            self.mw_control_seq([0, 0])
        elif mw_control == 'on':
            pass
        else:
            raise ValueError('unsupported mw_control parameter')

        self._start_device()
        time.sleep(self.time_pad)  # 先让激光和微波开几秒
        time.sleep(self.asg_dwell)
        counts = self.counter.getData().ravel().tolist()

        self.stop()

        # recover the asg control sequence for MW to be 'on'
        if mw_control == 'off':
            self.mw_control_seq(mw_seq_on)

        return counts


class PulseScheduler(FrequencyDomainScheduler):
    """
    Pulse-based ODMR manipulation scheduler
    """

    def __init__(self, *args, **kwargs):
        super(PulseScheduler, self).__init__(*args, **kwargs)
        self.name = 'Pulse ODMR Scheduler'

    def configure_odmr_seq(self, t_init, t_mw, t_read_sig, t_read_ref=None, inter_init_mw=3000, pre_read=200,
                           inter_mw_read=500, inter_readout=200, inter_period=200, N: int = 1000):
        """
        Wave form for single period:
            asg laser channel:
            -----                 -------------
            |   |                 |           |
            |   |-----------------|           |
            asg microwave channel:
                    -------------
                    |           |
            --------|           |--------------
            asg tagger acquisition channel:
                                   ---   ---
                                   | |   | |
            -----------------------| |---| |----
        All units for the parameters is 'ns'
        :param t_init: time span for laser initialization, e.g. 5000
        :param t_mw: time span for microwave actual operation in a ASG period, e.g. 800
        :param t_read_sig: time span for fluorescence signal readout, e.g. 400
        :param t_read_ref: time span for reference signal readout, e.g. 400
                            if the parameter is not assigned, means single-pulse readout
        :param inter_init_mw: time interval between laser initialization and MW operation pulses, e.g. 3000
        :param pre_read: previous time interval before Tagger readout and after laser readout pulses
        :param inter_mw_read: time interval between MW operation and laser readout pulses
        :param inter_readout: interval between single readout pulse and reference signal readout, e.g. 200
                            when t_read_ref is assigned, this parameter will play its role
        :param inter_period: interval between two neighbor periods, e.g. 200
        :param N: number of ASG operation periods
        """
        # unit: ns
        # total time for 'N' period, also for MW operation time at each frequency point
        if t_read_ref is not None:
            self.two_pulse_readout = True
            laser_seq = [t_init, inter_init_mw + t_mw + inter_mw_read,
                         pre_read + t_read_sig + inter_readout + t_read_ref + inter_period, 0]
            mw_seq = [0, t_init + inter_init_mw, t_mw,
                      inter_mw_read + pre_read + t_read_sig + inter_readout + t_read_ref + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw + inter_mw_read + pre_read, t_read_sig, inter_readout,
                          t_read_ref, inter_period]
            # apd_seq = [sum(tagger_seq[:-4]), sum(tagger_seq[-4:])]
        else:
            # single-pulse readout
            laser_seq = [t_init, inter_init_mw + t_mw + inter_mw_read, pre_read + t_read_sig + inter_period, 0]
            mw_seq = [0, t_init + inter_init_mw, t_mw, inter_mw_read + pre_read + t_read_sig + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw + inter_mw_read + pre_read, t_read_sig, inter_period]
            # apd_seq = [sum(tagger_seq[:-2], sum(tagger_seq[-2:]))]

        if self.mw_ttl == 0:
            #     # low-level effective
            mw_seq = utils.flip_sequence(mw_seq)

        self.download_asg_sequences(laser_seq, mw_seq, tagger_seq, N)

    def run_single_step(self, freq, power=None, mw_control='on') -> List[float]:
        """
        Single-frequency & single-power setting for running the scheduler
        :param freq: MW frequency, unit: Hz
        :param power: MW power, unit: dBm
        :param mw_control: 'on' or 'off'
        :return: 1-D array: [N,]
        """
        print('running for freq = {:.3f} GHz, power = {:.2f} dBm ...'.format(freq / C.giga, power))

        # set MW parameters
        self.configure_mw_paras(power, freq)

        # start sequence for time: N*t
        mw_seq_on = self.mw_control_seq()
        if mw_control == 'off':
            self.mw_control_seq([0, 0])
        elif mw_control == 'on':
            pass
        else:
            raise ValueError('unsupported mw_control parameter')

        self._start_device()
        time.sleep(0.5)  # 先让激光和微波开几秒
        time.sleep(self.asg_dwell)
        counts = self.counter.getData().ravel().tolist()
        self.stop()

        # recover the asg control sequence for MW to be 'on'
        if mw_control == 'off':
            self.mw_control_seq(mw_seq_on)

        return counts
