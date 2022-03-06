"""
Dynamical Decoupling detection (time-domain measurement)
---
1. configure pi/2 & pi pulses
2. scanning spin evolution time
"""

from odmactor.scheduler.base import Scheduler
import datetime
import time
import os
import numpy as np
import scipy.constants as C
from odmactor.utils import cal_contrast
import TimeTagger as tt
import pickle


class DynamicalDecouplingScheduler(Scheduler):
    def __init__(self, *args, **kwargs):
        super(DynamicalDecouplingScheduler, self).__init__(*args, **kwargs)
        self.name = 'Dynamical Decoupling Scheduler'

    def configure_odmr_seq(self, *args, **kwargs):
        pass
