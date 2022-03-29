from typing import List, Union
from odmactor.utils.asg import ASG8005

"""
ASG sequences example:
------------
asg_data1=[
    [100,700],
    [0,100,100,600],
    [0,200,100,500],
    [0,300,100,400],
    [0,400,100,300],
    [0,500,100,200],
    [0,600,100,100],
    [0,700,100,0]
]
"""


class ASG(ASG8005):
    def __init__(self):
        super(ASG, self).__init__()
        self.connect()

    def load_data(self, asg_data: List[List[Union[float, int]]]):
        """
        Connect ASG and download designed sequences data into it
        :param asg_data: ASG sequences for different channels
        """
        is_connected = super(ASG, self).connect()
        if is_connected == 1:
            super(ASG, self).download_ASG_pulse_data(asg_data, [len(row) for row in asg_data])
        else:
            raise ConnectionError('ASG not connected')

    def check_data(self, asg_data: List[List[Union[float, int]]]):
        return super(ASG, self).checkdata(asg_data, [len(row) for row in asg_data])

    def connect(self):
        super(ASG, self).connect()

    def start(self, count=1):
        super(ASG, self).start()

    def stop(self):
        super(ASG, self).stop()

    def close(self):
        super(ASG, self).close_device()
