from odmactor.scheduler import CWScheduler
import scipy.constants as C
import unittest
import numpy as np
from odmactor.utils.plotting import plot_freq_contrast
import  time

t_ns = 1e4
N = int(1e3)


# duration of single-frequency measurement: 10 ms

class SingleFrequencyTest(unittest.TestCase):
    def test_single_freq(self):
        f = 2.87 * C.giga
        p = 10

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

    def test_scanning_by_single(self):
        freq_start = 2.85 * C.giga
        freq_end = 2.89 * C.giga
        freq_step = 1 * C.mega
        p = 20

        scheduler = CWScheduler()
        scheduler.configure_mw_paras(p)
        scheduler.configure_odmr_seq(t_ns, N)
        scheduler.configure_tagger_counting()
        freqs = np.linspace(freq_start, freq_end, int((freq_end - freq_start) / freq_step + 1))
        counts = []
        for f in freqs:
            data_on = scheduler.run_single_step(p, f).mean()
            data_off = scheduler.run_single_step(p, f, mw_control=True).mean()
            counts.append([data_on, data_off])
        counts = np.array(counts)
        contrasts = np.abs(counts[:, 0] - counts[:, 1]) / counts[:, 0]

        print('Average counting:')
        print('MW on:')
        print(counts[:, 0].round(2))
        print('MW off:')
        print(counts[:, 1].round(2))
        print('Contrast:')
        print(counts.round(4))

        plot_freq_contrast(freqs, contrasts, fname='contrast-by-single')


if __name__ == '__main__':
    unittest.main()
