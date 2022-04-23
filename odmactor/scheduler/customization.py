"""
Customized Scheduler
"""
from odmactor.scheduler.base import Scheduler


class CustomizedScheduler(Scheduler):
    """
    Customized scheduler class
    """
    def __init__(self, *args, **kwargs):
        super(CustomizedScheduler, self).__init__(*args, **kwargs)
        self.name = 'Customized Scheduler'
