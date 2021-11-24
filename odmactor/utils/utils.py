import numpy as np


def cut_edge_zeros(arr):
    arr = arr[np.count_nonzero(arr, axis=1) != 0]  # strip rows
    arr = arr[:, np.count_nonzero(arr, axis=0) != 0]  # strip cols
    return arr


def cal_contrast(ls):
    """
    :param ls: size of [N*2,]
    :return: size of scalar
    """
    # Sequence location problems
    ls = np.array(ls)
    # 以下分母 max 是为了避免 CounterBetweenMarker Measurement 采集的数据错位
    contrast = np.abs(ls[::2].sum() - ls[1::2].sum()) / max(ls[1::2].sum(), ls[::2].sum())
    return contrast
    # if contrast.mean() > 1:
    #     contrast = 1 / contrast
    # return contrast.tolist()
