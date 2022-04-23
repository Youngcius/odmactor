"""
Spin manipulation classes
---

"""

from odmactor.scheduler.base import Scheduler


class SpinControlScheduler(Scheduler):
    """
    Spin control abstract class
    """

    def __init__(self, *args, **kwargs):
        super(SpinControlScheduler, self).__init__(*args, **kwargs)
        self.name = 'Base Spin Control Scheduler'
