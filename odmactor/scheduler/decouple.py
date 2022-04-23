"""
Dynamical Decoupling detection (time-domain measurement)
---
1. configure pi/2 & pi pulses
2. scanning spin evolution time
"""

from odmactor.scheduler.base import Scheduler


class DynamicalDecouplingScheduler(Scheduler):
    def __init__(self, *args, **kwargs):
        super(DynamicalDecouplingScheduler, self).__init__(*args, **kwargs)
        self.name = 'Dynamical Decoupling Scheduler'

    def configure_odmr_seq(self, *args, **kwargs):
        pass
