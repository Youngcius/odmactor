import numpy as np
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
                # 1) convert [0,0,...,0] into [0,0]
                if sum(seq) == 0:
                    asg_data[i] = [0, 0]
                else:
                    # 2) convert [a,b,c,d,0,0,...,0] into [a,b,c,d]
                    n = len(seq)
                    for j in range(n - 1, 0, -1):
                        if seq[j] > 0:
                            break
                        else:
                            asg_data[i].pop(j)
                    if len(asg_data[i]) % 2 != 0:
                        asg_data[i].append(0)
                    # 3) convert [a,b,c,0,d,e,f,0,g,h,i] into [a,b,c+d,e,f+g,h,i]
                    seq_mid = asg_data[i][1:-1]
                    seq_old = deepcopy(asg_data[i])
                    seq_new = []
                    n = len(asg_data[i])
                    if len(seq_mid) - np.count_nonzero(seq_mid) > 0:
                        seq_new.append(seq_old[0])
                        j = 1
                        while j < n - 1:
                            if seq_old[j] == 0:
                                seq_new[-1] += seq_old[j + 1]
                                j += 2
                            else:
                                seq_new.append(seq_old[j])
                                j += 1
                        seq_new.append(seq_old[-1])
                        if len(seq_new) % 2 != 0:
                            seq_new.append(0)
                        asg_data[i] = seq_new

        return asg_data

    def load_data(self, asg_data: List[List[int]]):
        """
        Connect ASG and download designed sequences data into it
        :param asg_data: ASG sequences for different channels
        """
        is_connected = super(ASG, self).connect()
        asg_data = self.normalize_data(asg_data)
        if is_connected == 1:
            return super(ASG, self).download_ASG_pulse_data(asg_data, [len(row) for row in asg_data])
        else:
            raise ConnectionError('ASG not connected')

    def check_data(self, asg_data: List[List[int]]):
        return super(ASG, self).checkdata(asg_data, [len(row) for row in asg_data])

    def connect(self):
        return super(ASG, self).connect()

    def start(self, count=1):
        return super(ASG, self).start()

    def stop(self):
        return super(ASG, self).stop()

    def close(self):
        return super(ASG, self).close_device()
