"""
Utils functions processing ASG sequences
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
from typing import List, Union
from operator import add
from copy import deepcopy
from matplotlib.figure import Figure


class SequenceString:
    """
    A OOP encapsulation: String representation of a sequence of pulses
    """

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


def sequences_to_figure(sequences: List[List[int]]) -> Figure:
    """
    Convert sequences (list of list) into a Figure instance
    """
    sequences = expand_to_same_length(sequences)
    N = len(sequences)  # num_channels
    idx_exist = [i for i, l in enumerate(sequences) if sum(l) > 0]
    n = len(idx_exist)  # effective number of channels
    channels = ['ch {}'.format(i + 1) for i in range(N)]
    seq_eff = [sequences[i] for i in idx_exist]
    gcd = reduce(math.gcd, list(map(int, reduce(add, seq_eff))))
    for i in range(n):
        seq_eff[i] = [int(t / gcd) for t in seq_eff[i]]
    baselines = []
    levels = []

    seq_all = [[] for i in range(N)]
    length = sum(seq_eff[0])
    for i in range(N):
        if i in idx_exist:
            seq_all[i] = seq_eff[idx_exist.index(i)]
        else:
            seq_all[i] = [0, length]  # 0 个 '1', length 个 '0'

    for i in range(N):
        # 0,1,2,3,...,N-1
        level = []
        j = 0
        l = len(seq_all[i])
        while j < l:
            level += [1] * seq_all[i][j] + [0] * seq_all[i][j + 1]
            j += 2

        b = 1.2 * i
        level = [lev + b for lev in level]
        baselines.append(b)
        levels.append(level)

    fig = plt.figure(figsize=(14, 2 * len(idx_exist)))
    for i, ch in enumerate(channels):
        plt.stairs(levels[i], baseline=baselines[i] - 0.03, label=ch, fill=True)
    plt.title('Sequences', fontsize=20)
    plt.xlabel('time ({} ns)'.format(int(gcd)), fontsize=15)
    plt.yticks(baselines, channels, fontsize=13)

    plt.xlim(0, max([sum(s) for s in seq_eff]))
    plt.ylim(-0.1, max(levels[-1]) + 0.1)
    plt.xticks(fontsize=13)
    return fig


# def seq_to_fig(seq: List[List[Union[float, int]]]) -> Figure:
#     """
#     Convert sequences (list of list) into a Figure instance
#     """
#     idx_exist = [i for i, l in enumerate(seq) if sum(l) > 0]
#     n = len(idx_exist)  # num_channels
#     channels = ['ch {}'.format(i + 1) for i in idx_exist]
#     seq_eff = [seq[i] for i in idx_exist]
#     gcd = reduce(math.gcd, list(map(int, reduce(concat, seq_eff))))
#     for i in range(n):
#         seq_eff[i] = [int(t / gcd) for t in seq_eff[i]]
#     baselines = []
#     levels = []

#     for i in range(n):
#         # 0,1,2,3,...
#         level = []
#         j = 0
#         l = len(seq_eff[i])
#         while j < l:
#             level += [1] * seq_eff[i][j] + [0] * seq_eff[i][j + 1]
#             j += 2

#         b = 1.1 * i
#         level = [lev + b for lev in level]
#         baselines.append(b)
#         levels.append(level)
#     fig = plt.figure(figsize=(14, 2 * len(idx_exist)))
#     for i, ch in enumerate(channels):
#         plt.stairs(levels[i], baseline=baselines[i], label=ch)
#     plt.legend(loc='upper left')
#     plt.title('Sequences')
#     plt.ylabel('channel')
#     plt.xlabel('time ({} ns)'.format(int(gcd)))
#     plt.xlim(0, max([sum(s) for s in seq_eff]))
#     plt.ylim(-0.1, max(levels[-1]) + 0.1)
#     plt.yticks([])
#     return fig


def sequences_to_string(sequences: List[List[int]]) -> str:
    """
    Convert sequences (list of list) into a string
    """
    sequences = expand_to_same_length(sequences)
    idx_exist = [i for i, l in enumerate(sequences) if sum(l) > 0]

    # flatten; float --> integer; calculate gcd
    gcd = reduce(math.gcd, list(map(int, reduce(add, sequences))))

    str_dict = {'channel {}'.format(i + 1): SequenceString() for i in idx_exist}

    for i in idx_exist:
        length = len(sequences[i])
        j = 0
        while j < length:
            # pair of a high pulse and a low pulse
            high_width = int(sequences[i][j] / gcd)
            low_width = int(sequences[i][j + 1] / gcd)
            str_dict['channel {}'.format(i + 1)].append_high_pulse(high_width)
            str_dict['channel {}'.format(i + 1)].append_low_pulse(low_width)
            j += 2

    str_list = ['\n'.join([k, str(v)]) for k, v in str_dict.items()]
    return '\n\n'.join(str_list)


def expand_to_same_length(sequences: List[List[int]]) -> List[List[int]]:
    """
    Expand eac sequence to the same length, by calculating the LCM of all sequences lengths
    """
    sequences_expanded = deepcopy(sequences)
    lengths = [int(sum(seq)) for seq in sequences]
    if len(np.unique(lengths)) == 1:
        return sequences_expanded
    else:
        tm = np.lcm.reduce(lengths)
        for i, t in enumerate(lengths):
            sequences_expanded[i] *= int(tm / t)
        return sequences_expanded


def flip_sequence(seq: list) -> list:
    """
    Flip the control sequence
    i.e., high-level effective <---> low level effective
    """
    if seq[0] == 0:
        if seq[-1] == 0:
            return seq[1:-1]
        else:
            return seq[1:] + [0]
    else:
        if seq[-1] == 0:
            return [0] + seq[:-1]
        else:
            return [0] + seq + [0]
