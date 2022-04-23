import time
from odmactor.utils import sequence
seq = [
    [3000.0, 3000, 2000, 0],
    [0, 4000.0, 1000, 3000],
    [0, 0],
    [1000000, 1000000],
    [0, 6000.0, 700, 1300],
    [0, 0],
    [0, 0],
    [1000000, 1000000]
]

seq_same = sequence.expand_to_same_length(seq)

import matplotlib.pyplot as plt
sequence.sequences_to_figure(seq_same)
plt.show()






# ===============================================
#
#
#
#
# seq1 = [
#     [3000.0, 3000, 900, 0],
#     [0, 4000.0, 1000, 1900],
#     [0, 0],
#     [1000, 0],
#     [0, 6000.0, 700, 200],
#     [0, 0],
#     [0, 0],
#     [1000, 0]
# ]
#
# seq0 = [
#     [3000.0, 3000, 900, 0],
#     [0, 4000.0, 1000, 1900],
#     [0, 0],
#     [0, 0],
#     [0, 6000.0, 700, 200],
#     [0, 0],
#     [0, 0],
#     [0, 0]
# ]
#
# from odmactor.instrument import ASG
#
# asg = ASG()
# print(asg.connect(), asg.load_data(seq), asg.load_data(seq0), asg.load_data(seq1))
# asg.load_data(
#     [
#         [1000, 1000],
#         [1000, 1000],
#         [1000, 1000],
#         [1000, 1000],
#         [1000, 1000],
#         [1000, 1000],
#         [1000, 1000],
#         [1000, 1000],
#
#     ]
# )
#
# print(asg.normalize_data(seq))
# print(asg.normalize_data(seq0))
# print(asg.normalize_data(seq1))
#
# time.sleep(100)
