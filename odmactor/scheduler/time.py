"""
Time-domain ODMR detecting classes
---
1. Ramsey detecction
    1) configure pi/2 pulse
    2) scanning free precession time
2. Rabi oscillation
    1) fix the power of MW & freq
    2) set the range of MW duration time, for scanning
    3) repeat execution at each duration time point
3. T1 Relaxometry
    1) configure pi pulse
    2) scanning free precession time
"""

from odmactor.scheduler.base import TimeDomainScheduler
import time
import scipy.constants as C
from typing import List
from odmactor import utils


class RamseyScheduler(TimeDomainScheduler):
    """
    Ramsey detection scheduler
    """

    def __init__(self, *args, **kwargs):
        super(RamseyScheduler, self).__init__(*args, **kwargs)
        self.name = 'Ramsey Scheduler'

    def _gene_detect_seq(self, t_free):
        """
        Generate Ramsey sequences and download it to ASG
        :param t_free: free precession time (time duration between two MW pulse)
        """
        t_init, t_mw = self._cache['t_init'], self._cache['t_mw']
        inter_init_mw, inter_mw_read = self._cache['inter_init_mw'], self._cache['inter_mw_read']
        t_read_sig, t_read_ref = self._cache['t_read_sig'], self._cache['t_read_ref']
        inter_readout, pre_read = self._cache['inter_readout'], self._cache['pre_read']
        inter_period = self._cache['inter_period']
        N = self._cache['N']

        # generate ASG wave forms
        if self.two_pulse_readout:
            laser_seq = [t_init, inter_init_mw + t_mw * 2 + t_free + inter_mw_read,
                         pre_read + t_read_sig + inter_readout + t_read_ref + inter_period, 0]
            mw_seq = [0, t_init + inter_init_mw, t_mw, t_free, t_mw,
                      inter_mw_read + pre_read + t_read_sig + inter_readout + t_read_ref + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw * 2 + t_free + inter_mw_read + pre_read, t_read_sig,
                          inter_readout, t_read_ref, inter_period]
        else:
            # single-pulse readout (without reference readout)
            laser_seq = [t_init, inter_init_mw + t_mw * 2 + t_free + inter_mw_read,
                         pre_read + t_read_sig + inter_period, 0]
            mw_seq = [0, t_init + inter_init_mw, t_mw, t_free, t_mw,
                      inter_mw_read + pre_read + t_read_sig + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw * 2 + t_free + inter_mw_read + pre_read, t_read_sig,
                          inter_period]

        if self.mw_ttl == 0:
            mw_seq = utils.flip_sequence(mw_seq)

        self.download_asg_sequences(laser_seq, mw_seq, tagger_seq, N)

    def configure_odmr_seq(self, t_init, t_read_sig, t_read_ref=None, inter_init_mw=1000, inter_mw_read=200,
                           inter_readout=200,
                           pre_read=50, inter_period=200, N: int = 1000):
        """
        Wave form for single period:
            laser (no asg control sequence):
            -----                       ---------
            |   |                       |       |
            |   |-----------------------|       |
            asg microwave channel:
                  ---               ---
                  | |  (<------->)  | |
            ------| |---------------| |----------
            asg tagger acquisition channel:
                                        ---   ---
                                        | |   | |
            ----------------------------| |---| |
        All units for the parameters is 'ns'
        :param t_init: time for laser initialization
        :param t_read_sig: time span for fluorescence signal readout
        :param t_read_ref: time span for reference signal readout
        :param inter_init_mw: time interval between laser initialization and MW operation pulses, e.g. 1000
        :param inter_readout: time span for interval
        :param inter_period: interval between two neighbor periods, e.g. 200
        :param N: number of ASG operation periods for each detection point
        """
        t_mw = self.pi_pulse['time'] / 2 / C.nano  # pi/2 pulse time duration, s --> ns
        self._cache = {
            't_init': t_init,
            't_mw': t_mw,
            't_read_sig': t_read_sig,
            't_read_ref': t_read_ref,
            'inter_init_mw': inter_init_mw,
            'inter_mw_read': inter_mw_read,
            'inter_readout': inter_readout,
            'pre_read': pre_read,
            'inter_period': inter_period,
            'N': N,
        }
        self._asg_conf['N'] = N
        if t_read_ref is not None:
            self.two_pulse_readout = True

    def run_single_step(self, t_free, power=None) -> List[float]:
        """
        Single-interval & single-power setting for running the scheduler
        :param power: MW power, unit: dBm
        :param t_free: free precession time, unit: ns
        :return: 1-D array: [N,]
        """
        print('running for time interval = {:.3f} ns, power = {:.2f} dBm ...'.format(t_free, power))

        # set MW parameters
        self.configure_mw_paras(power, regulate_pi=True)

        # generate ASG sequences
        self._gene_detect_seq(t_free)

        # start sequence for time: N*t
        self._start_device()
        time.sleep(0.5)  # 先让激光和微波开几秒
        time.sleep(self.asg_dwell)
        counts = self.counter.getData().ravel().tolist()
        self.stop()

        return counts


class RabiScheduler(TimeDomainScheduler):
    """
    Ramsey detection scheduler
    """

    def __init__(self, *args, **kwargs):
        super(RabiScheduler, self).__init__(*args, **kwargs)
        self.name = 'Rabi Scheduler'

    def _gene_detect_seq(self, t_mw):
        """
        Generate Rabi detection sequences and download it to ASG
        :param t_mw: free precession time (time duration between two MW pulse)
        """
        t_init, inter_init_mw = self._cache['t_init'], self._cache['inter_init_mw']
        t_read_sig, t_read_ref = self._cache['t_read_sig'], self._cache['t_read_ref']
        inter_readout, inter_period = self._cache['inter_readout'], self._cache['inter_period']
        inter_mw_read, pre_read = self._cache['inter_mw_read'], self._cache['pre_read']
        N = self._cache['N']

        # generate ASG wave forms
        if self.two_pulse_readout:
            laser_seq = [t_init, inter_init_mw + t_mw + inter_mw_read,
                         pre_read + t_read_sig + inter_readout + t_read_ref + inter_period, 0]
            mw_seq = [0, t_init + inter_init_mw, t_mw,
                      inter_mw_read + pre_read + t_read_sig + inter_readout + t_read_ref + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw + inter_mw_read + pre_read, t_read_sig, inter_readout,
                          t_read_ref, inter_period]
        else:
            # single-pulse readout (without reference readout)
            laser_seq = [t_init, inter_init_mw + t_mw + inter_mw_read, pre_read + t_read_sig + inter_period, 0]
            mw_seq = [0, t_init + inter_init_mw, t_mw, inter_mw_read + pre_read + t_read_sig + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw + inter_mw_read + pre_read, t_read_sig, inter_period]

        if self.mw_ttl == 0:
            mw_seq = utils.flip_sequence(mw_seq)

        self.download_asg_sequences(laser_seq, mw_seq, tagger_seq, N)

    def configure_odmr_seq(self, t_init, t_read_sig, t_read_ref=None, inter_init_mw=1000, inter_mw_read=100,
                           pre_read=200, inter_readout=200, inter_period=200, N: int = 1000):
        """
        Wave form for single period:
            laser (no asg control sequence):
            -----                        ---------
            |   |                        |       |
            |   |------------------------|       |
            asg microwave channel (variable duration):
                  ---------------------
                  |  (<----------->)  |
            ------|                   |----------
            asg tagger acquisition channel:
                                         ---   ---
                                         | |   | |
            -----------------------------| |---| |
        All units for the parameters is 'ns'
        :param t_init: time for laser initialization
        :param t_read_sig: time span for fluorescence signal readout
        :param t_read_ref: time span for reference signal readout
        :param inter_init_mw: time interval between laser initialization and MW operation pulses, e.g. 1000
        :param inter_mw_read: time interval between MW operation and readout laser pulses
        :param pre_read: previous time interval before Tagger readout and after laser readout pulses
        :param inter_readout: time span for interval
        :param inter_period: interval between two neighbor periods, e.g. 200
        :param N: number of ASG operation periods for each detection point
        """
        self._cache = {
            't_init': t_init,
            't_read_sig': t_read_sig,
            't_read_ref': t_read_ref,
            'inter_init_mw': inter_init_mw,
            'inter_mw_read': inter_mw_read,
            'pre_read': pre_read,
            'inter_readout': inter_readout,
            'inter_period': inter_period,
            'N': N,
        }
        self._asg_conf['N'] = N
        if t_read_ref is not None:
            self.two_pulse_readout = True

    def run_single_step(self, t_mw, power=None) -> List[float]:
        """
        Single-interval & single-power setting for running the scheduler
        :param power: MW power, unit: dBm
        :param t_mw: free precession time, unit: ns
        :return: 1-D array: [N,]
        """
        print('running for MW operation time = {:.3f} ns, power = {:.2f} dBm ...'.format(t_mw, power))

        # set MW parameters
        self.configure_mw_paras(power)

        # generate ASG sequences
        self._gene_detect_seq(t_mw)

        # start sequence for time: N*t
        self._start_device()
        time.sleep(0.5)  # 先让激光和微波开几秒
        time.sleep(self.asg_dwell)
        counts = self.counter.getData().ravel().tolist()
        self.stop()

        return counts


class RelaxationScheduler(TimeDomainScheduler):
    """
    T1 Relaxation measurement Scheduler
    """

    def __init__(self, *args, **kwargs):
        super(RelaxationScheduler, self).__init__(*args, **kwargs)
        self.name = 'T1 Relaxation Scheduler'
        if 'ms' in kwargs.keys():
            self.ms = 0
        else:
            self.ms = 1

    def _gene_detect_seq(self, t_free):
        """
        Generate T1 Relaxation detection sequences and download it to ASG
        :param t_free: free precession time (time duration of after MW pi pulse)
        """
        t_init, inter_init_mw, t_mw = self._cache['t_init'], self._cache['inter_init_mw'], self._cache['t_mw']
        t_read_sig, t_read_ref = self._cache['t_read_sig'], self._cache['t_read_ref']
        inter_readout, inter_period = self._cache['inter_readout'], self._cache['inter_period']
        N = self._cache['N']

        if self.ms == 1:
            if self.two_pulse_readout:
                # meaning t_read_ref is not None:
                laser_seq = [t_init, inter_init_mw + t_mw + t_free, t_read_sig + inter_readout + t_read_ref,
                             inter_period]
                mw_seq = [0, t_init + inter_init_mw, t_mw,
                          t_free + t_read_sig + inter_readout + t_read_ref + inter_period]
                tagger_seq = [0, t_init + inter_init_mw + t_mw + t_free, t_read_sig, inter_readout, t_read_ref,
                              inter_period]
            else:
                # single-pulse readout (without reference readout)
                laser_seq = [t_init, inter_init_mw + t_mw + t_free, t_read_sig, inter_period]
                mw_seq = [0, t_init + inter_init_mw, t_mw, t_free + t_read_sig + inter_period]
                tagger_seq = [0, t_init + inter_init_mw + t_mw + t_free, t_read_sig, inter_period]
        else:
            print('ms == 0')
            if self.two_pulse_readout:
                laser_seq = [t_init, t_free, t_read_sig + inter_readout + t_read_ref, inter_period]
                mw_seq = [0, sum(laser_seq)]
                tagger_seq = [0, t_init + t_free, t_read_sig, inter_readout, t_read_ref, inter_period]
            else:
                laser_seq = [t_init, t_free, t_read_sig, inter_period]
                mw_seq = [0, sum(laser_seq)]
                tagger_seq = [0, t_init + t_free, t_read_sig, inter_period]

        if self.mw_ttl == 0:
            mw_seq = utils.flip_sequence(mw_seq)

        self.download_asg_sequences(laser_seq, mw_seq, tagger_seq, N)

    def configure_odmr_seq(self, t_init, t_read_sig, t_read_ref=None, inter_init_mw=10000, inter_readout=200,
                           inter_period=200, N: int = 10000):
        """
        Wave form for single period:
            laser (no asg control sequence):
            -----                      ---------
            |   |                      |       |
            |   |----------------------|       |
            asg microwave channel (variable duration):
                  -----
                  |   | (<--------->)    ms = 1
            ------|   |------------------------
            or
                 (<----------------->)  ms = 0
            -----------------------------------
            asg tagger acquisition channel:
                                       ---   ---
                                       | |   | |
            ---------------------------| |---| |
        All units for the parameters is 'ns'
        :param t_init: time for laser initialization
        :param t_read_sig: time span for fluorescence signal readout
        :param t_read_ref: time span for reference signal readout
        :param inter_init_mw: time interval between laser initialization and MW operation pulses, e.g. 1000
        :param inter_readout: time span for interval
        :param inter_period: interval between two neighbor periods, e.g. 200
        :param N: number of ASG operation periods for each detection point
        """
        t_mw = self.pi_pulse['time'] / C.nano  # pi/2 pulse time duration, s --> ns
        self._cache = {
            't_init': t_init,
            't_mw': t_mw,
            't_read_sig': t_read_sig,
            't_read_ref': t_read_ref,
            'inter_init_mw': inter_init_mw,
            'inter_readout': inter_readout,
            'inter_period': inter_period,
            'N': N,
        }
        self._asg_conf['N'] = N
        if t_read_ref is not None:
            self.two_pulse_readout = True

    def run_single_step(self, t_free, power=None) -> List[float]:
        """
        Single-interval & single-power setting for running the scheduler
        :param power: MW power, unit: dBm
        :param t_free: free precession time (time duration of after MW pi pulse)
        :return: 1-D array: [N,]
        """
        print('running for relaxation time = {:.3f} ns, power = {:.2f} dBm ...'.format(t_free, power))

        # set MW parameters
        self.configure_mw_paras(power, regulate_pi=True)

        # generate ASG sequences
        self._gene_detect_seq(t_free)

        # start sequence for time: N*t
        self._start_device()
        time.sleep(0.5)  # 先让激光和微波开几秒
        time.sleep(self.asg_dwell)
        counts = self.counter.getData().ravel().tolist()
        self.stop()

        return counts
