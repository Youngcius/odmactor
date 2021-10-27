from odmactor.scheduler.base import Scheduler
import datetime
import time
import os
import numpy as np
import scipy.constants as C
from odmactor.utils import cal_contrast
import TimeTagger as tt
import pickle

"""
Rabi oscillation measurement
---
1. fix the power of MW & freq
2. set the range of MW duration time, for sanning
3. repeat execution at each duration time point
"""


class RabiScheduler(Scheduler):
    """
    Ramsey detection scheduler
    """

    def __init__(self, *args, **kwargs):
        super(RabiScheduler, self).__init__()
        self.name = 'Rabi ODMR Scheduler'
        self._times = np.array([])

    def configure_mw_time_range(self, tmin, tmax, tstep=None, tnum=None):
        """
        Configure MW scanning time range for each time duration point, unit: ns
        :param tmin: minimum time duration
        :param tmax: maximum time duration
        :param tstep: time duration step
        :param tnum: number of time duration points
        """
        if tstep is not None:
            self._times = np.arange(tmin, tmax + 0.1, tstep)
        elif tnum is not None:
            self._times = np.linspace(tmin, tmax, tnum)
        else:
            raise TypeError('Please input an efficient "tstep" or "tnum" parameter')

    def configure_part_odmr_seq(self, t_init, t_read_sig, t_read_ref, N: int = 100, t_interval=20):
        """
        :param t_init: time for laser initialization
        :param t_read_sig: time span for fluorescence signal readout
        :param t_read_ref: time span for reference signal readout
        :param N: number of ASG operation periods
        :param t_interval: time span for interval
        """
        self._cache = {
            't_init': t_init,
            't_read_sig': t_read_sig,
            't_read_ref': t_read_ref,
            'N': N,
            't_interval': t_interval,
        }

    def configure_odmr_seq(self, t_init, t_mw, t_read_sig, t_read_ref, t_interval=20):
        """
        Wave form for single period:
            laser (no asg control sequence):
            -----                   ---------
            |   |                   |       |
            |   |-------------------|       |
            asg microwave channel (variable duration):
                ---------------------
                |                   |
            ----|                   |--------
            asg tagger acquisition channel:
                                    ---   ---
                                    | |   | |
            ------------------------| |---| |
        All units for the parameters is 'ns'
        :param t_init: time for laser initialization
        :param t_mw: time for microwave operation
        :param t_read_sig: time span for fluorescence signal readout
        :param t_read_ref: time span for reference signal readout
        :param t_interval: time span for interval
        """
        # unit: ns
        # self._asg_conf['t'] = t * C.nano

        # generate ASG wave forms
        laser_seq = [t_init, 2 * t_interval + t_mw, t_read_sig + t_interval + t_read_ref, 0]
        mw_seq = [0, t_init + t_interval, t_mw, t_read_sig + t_read_ref + 2 * t_interval]
        tagger_seq = [0, t_init + t_mw + 2 * t_interval, t_read_sig, t_interval, t_read_ref, 0]

        idx_laser_channel = self.channel['laser'] - 1
        idx_mw_channel = self.channel['w'] - 1
        idx_tagger_channel = self.channel['tagger'] - 1
        idx_apd_channel = self.channel['apd'] - 1

        self._asg_sequences[idx_laser_channel] = laser_seq
        self._asg_sequences[idx_mw_channel] = mw_seq
        self._asg_sequences[idx_tagger_channel] = tagger_seq
        self._asg_sequences[idx_apd_channel] = tagger_seq

        # connect & download pulse data
        self.asg_connect_and_download_data(self._asg_sequences)

    def _start_device(self, *args, **kwargs):
        # run MW firstly
        self._mw_instr.write_bool('OUTPUT:STATE', True)
        if self.mw_exec_mode == 'scan-center-span' or self.mw_exec_mode == 'scan-start-stop':
            self._mw_instr.write_str('SWE:FREQ:EXEC')  # trigger the sweep
        else:
            self._mw_instr.write_str()  # 在外在激励的作用下执行N个周期 TODO

        # run ASG then
        self._asg.start()

        # execute Measurement instance
        self.counter = tt.CountBetweenMarkers(self.tagger, self.tagger_input['apd'],
                                              begin_channel=self.tagger_input['asg'],
                                              end_channel=-self.tagger_input['asg'],
                                              n_values=self._asg_conf['N'] * 2)
        self.counter.startFor(int((self.time_total + 5) / C.pico))  # parameter unit: ps

    def _acquire_data(self, duration: float, *args, **kwargs):
        time.sleep(self.mw_dwell)
        data = self.counter.getData()  # [N*2,]
        return data.tolist()

    def run(self):
        """
        iteration: scanning MW operation time (ASG MW channel):
            1) start devices
            2) acquire data
        save data result
        """
        t_init, t_read_sig, t_read_ref = self._cache['t_int'], self._cache['t_read_sig'], self._cache['t_read_ref']
        t_interval, N = self._cache['t_interval'], self._cache['N']
        self._asg_conf['N'] = N

        data = []
        for t_mw in self._times:
            print('scanning MW time duration: {:.2f} ns ...'.format(t_mw))
            self._asg_conf['t'] = t_init + t_mw + 3 * t_interval + t_read_sig + t_read_ref
            self.mw_dwell = self._cache['N'] * (self._asg_conf['t']) * C.nano  # unit: s
            self.configure_odmr_seq(t_init, t_mw, t_read_sig, t_read_ref, t_interval)
            self._start_device()
            data.append(self._acquire_data(self.mw_dwell))
            self.stop()

        # calculate contrast & save result
        contrast_list = [cal_contrast(ls) for ls in data]
        contrast = [np.mean(item) for item in contrast_list]  # average contrast for each t_mw duration point
        self._result = np.array([self._times, contrast])
        self._result_detail = {
            'times': self._times,
            'contrast': contrast,
            'origin_data': data,
            'contrast_list': contrast_list
        }
        fname = os.path.join(self.output_dir,
                             'Rabi-ODMR-result-{}-{}'.format(datetime.date.today(), round(time.time() / 120)))
        np.savetxt(fname + '.txt', self._result)
        with open(fname + '.pkl', 'wb') as f:
            pickle.dump(self._result_detail, f)
        print('data has been saved into {}'.format(fname))
