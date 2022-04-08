from copy import deepcopy
from typing import List
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

    def normalize_data(self, sequences: List[List[int]]) -> List[List[int]]:
        """
        Normalize sequences to make it acceptable data for ASG loading
        """
        asg_data = deepcopy(sequences)
        if not self.check_data(asg_data):
            for i, seq in enumerate(asg_data):
                if sum(seq) == 0:
                    asg_data[i] = [0, 0]
                else:
                    n = len(seq)
                    for j in range(n - 1, 0, -1):
                        if seq[j] > 0:
                            break
                        else:
                            asg_data[i].pop(j)
                    if len(asg_data[i]) % 2 != 0:
                        asg_data[i].append(0)
        return asg_data

    def load_data(self, asg_data: List[List[int]]):
        """
        Connect ASG and download designed sequences data into it
        :param asg_data: ASG sequences for different channels
        """
        is_connected = super(ASG, self).connect()
        if is_connected == 1:
            super(ASG, self).download_ASG_pulse_data(asg_data, [len(row) for row in asg_data])
        else:
            raise ConnectionError('ASG not connected')

    def check_data(self, asg_data: List[List[int]]):
        return super(ASG, self).checkdata(asg_data, [len(row) for row in asg_data])

    def connect(self) -> bool:
        return super(ASG, self).connect()

    def start(self, count=1):
        return super(ASG, self).start()

    def stop(self):
        return super(ASG, self).stop()

    def close(self):
        return super(ASG, self).close_device()
