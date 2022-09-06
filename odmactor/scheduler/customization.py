"""
User-customized Scheduler
"""
import time
import scipy.constants as C
from typing import List
from .base import Scheduler


class CustomizedScheduler(Scheduler):
    """
    Customized scheduler class
    """

    def __init__(self, *args, **kwargs):
        super(CustomizedScheduler, self).__init__(*args, **kwargs)
        self.name = 'Customized Scheduler'

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
