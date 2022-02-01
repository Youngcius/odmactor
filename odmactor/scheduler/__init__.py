"""
Priority to be implemented: CW, Pulse, Ramsey, Rabi
TODO: Dynamical Decoupling, Hahn echo, T1 relaxation
Reference: https://arxiv.org/pdf/1910.00061.pdf
"""

from odmactor.scheduler.base import Scheduler
from odmactor.scheduler.odmr import CWScheduler, PulseScheduler
from odmactor.scheduler.ramsey import RamseyScheduler
from odmactor.scheduler.spin import SpinControlScheduler
