import math
from functools import reduce
from typing import List
from operator import concat
import numpy as np

"""
Utils functions processing ASG sequences
"""


class SequenceString:
    def __init__(self, height: int = 3, unit_width=1):
        self.strings = [''] * height
        self.unit_width = unit_width

    def append_low_pulse(self, width):
        l = len(self.strings)
        for i in range(l - 1):
            # for string in self.strings[:-1]:
            self.strings[i] += ' ' * width * self.unit_width
        self.strings[-1] += '-' * width * self.unit_width

    def append_high_pulse(self, width):
        l = len(self.strings)
        if width == 0:
            pass
        else:
            for i in range(1, l):
                # for string in self.strings[1:]:
                self.strings[i] += '|'
                self.strings[i] += ' ' * width * self.unit_width
                self.strings[i] += '|'
            self.strings[0] += '-' * width * self.unit_width + '-' * 2

    def __str__(self):
        return '\n'.join(self.strings)


def seq_to_str(seq: List[List[float]]):
    """
    Convert sequences (list of list) into a string
    """
    idx_exist = [i for i, l in enumerate(seq) if sum(l) > 0]

    # flatten, to integer, calculate gcd
    gcd = reduce(math.gcd, list(map(int, reduce(concat, seq))))

    str_dict = {'channel {}'.format(i + 1): SequenceString() for i in idx_exist}

    for i in idx_exist:
        length = len(seq[i])
        j = 0
        while j < length:
            # pair of a high pulse and a low pulse
            high_width = int(seq[i][j] / gcd)
            low_width = int(seq[i][j + 1] / gcd)
            str_dict['channel {}'.format(i + 1)].append_high_pulse(high_width)
            str_dict['channel {}'.format(i + 1)].append_low_pulse(low_width)
            j += 2

    str_list = ['\n'.join([k, str(v)]) for k, v in str_dict.items()]
    return '\n\n'.join(str_list)


if __name__ == '__main__':
    # 测试
    #             -----               ---------
    #             |   |               |       |
    #             |   |---------------|       |----
    #             asg microwave channel:
    #                     -------------
    #                     |           |
    #             --------|           |------------
    #             asg tagger acquisition channel:
    #                                 ---   ---
    #                                 | |   | |
    #             --------------------| |---| |----
    asgdata = [
        [100, 100 + 200, 100 + 20 + 100 + 30, 0],
        [0, 100 + 100, 200, 100 + 20 + 100 + 30],
        [0, 100 + 100 + 200, 100, 20, 100, 30],
        [0, 0],
        [0, 100 + 100 + 200, 100, 20, 100, 30],
        [0, 0],
        [0, 0],
        [0, 0],
    ]

    print(seq_to_str(asgdata))
