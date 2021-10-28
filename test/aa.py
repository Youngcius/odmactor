from odmactor.scheduler import CWScheduler
import scipy.constants as C
import unittest
import numpy as np
from odmactor.utils.plotting import plot_freq_contrast
import  time

t_ns = 1e4
N = int(1e3)

def single_freq():
    f = 2.87 * C.giga
    p = 20

    scheduler = CWScheduler()
    scheduler.configure_mw_paras(p)
    scheduler.configure_odmr_seq(t_ns, N)
    scheduler.configure_tagger_counting()  # default: Counter
    data_on = scheduler.run_single_step(p, f).mean()
    time.sleep(0.5)
    data_off = scheduler.run_single_step(p, f, mw_control=True).mean()

    # 如果是 4 Mps 的计数，10 us 均值 ~ 40
    print('Average counting:')
    print('MW on: {:.2f}'.format(data_on), 'MW off: {:.3f}'.format(data_off))
    print('contrast: {:.4f}'.format(np.abs(data_on - data_off) / data_on))

if __name__ == '__main__':
    single_freq()