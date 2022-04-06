"""
Abstract interfaces of instruments, i.e., Laser, Microwave, ASG, ...
"""

from .asg import ASG
from .laser import Laser
from .microwave import Microwave
from .lockin import LockInAmplifier