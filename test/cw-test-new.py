import matplotlib.pyplot as plt
from odmactor.scheduler import CWScheduler
import scipy.constants as C
import unittest
import numpy as np
from odmactor.utils.plotting import plot_freq_contrast
import time

t_ns = 1e5
N = int(1e4)


# duration of single-frequency measurement: 10 ms

class SingleFrequencyTest():
    def test_single_freq(self, p, f):
        f = f * C.giga
        p = p
        scheduler = CWScheduler()
        scheduler.configure_mw_paras(p)
        scheduler.configure_odmr_seq(t_ns, N)
        scheduler.configure_tagger_counting()  # default: Counter
        time.sleep(0.5)
        for i in range(3):
            data_on = scheduler.run_single_step(p, f).mean()
            data_off = scheduler.run_single_step(p, f, mw_control='off').mean()
            # print()
            # 如果是 4 Mps 的计数，1 us 均值 ~ 4
            # print('Average counting:')
            print('MW on: {:.4f}'.format(data_on), '   MW off: {:.4f}'.format(data_off))
            print('contrast: {:.4f}'.format(np.abs(data_on - data_off) / data_on))
        print()

    def test_scanning_by_single(self):
        freq_start = 2.85 * C.giga
        freq_end = 2.89 * C.giga
        freq_step = 1 * C.mega
        p = 10

        scheduler = CWScheduler()
        scheduler.configure_mw_paras(p)
        scheduler.configure_odmr_seq(t_ns, N)
        scheduler.configure_tagger_counting()  # default: Counter
        freqs = np.linspace(freq_start, freq_end, int((freq_end-freq_start)/freq_step + 1))
        for f in freqs:
            print('scanning freq {:.3f} GHz (by single)'.format(f/C.giga))

            data_on = scheduler.run_single_step(p, f).mean()
        time.sleep(0.5)
        for i in range(5):
            data_off = scheduler.run_single_step(p, f, mw_control='off').mean()
            print()
            # 如果是 4 Mps 的计数，1 us 均值 ~ 4
            # print('Average counting:')
            print('MW on: {:.4f}'.format(data_on), '   MW off: {:.4f}'.format(data_off))
            print('contrast: {:.4f}'.format(np.abs(data_on - data_off) / data_on))

    def test_scanning(self):
        freq_start = 2.85 * C.giga
        freq_end = 2.89 * C.giga
        freq_step = 1 * C.mega
        p = 10
        # 通道、微波、序列、counter
        scheduler = CWScheduler()
        scheduler.configure_mw_paras(p)
        scheduler.set_mw_freqs(freq_start, freq_end, freq_step)
        scheduler.configure_odmr_seq(t_ns, N)
        scheduler.configure_tagger_counting()

        scheduler.run_scanning('on')
        res_on = scheduler.result  # [freqs, counts]

        scheduler.run_scanning('off')
        res_off = scheduler.result
        scheduler.close()
        contrasts = [abs(c_off - c_on) / c_on for c_on, c_off in zip(res_on[1], res_off[1])]

        plt.figure()
        plt.plot(*res_on, label='MW on')
        plt.plot(*res_off, label='MW off')
        plt.legend()
        plt.title('Counts comparison')
        plt.savefig('counts_on_off (CW)', dpi=400)

        plot_freq_contrast(res_on[0], contrasts, fname='contrasts-CW')

        print('Average counting:')
        print('MW on:')
        print(np.round(res_on, 2))
        print('MW off:')
        print(np.round(res_off, 2))
        print('Contrast:')
        print(np.round(contrasts))
        #

        np.savetxt('counts-cw-new.txt', res_on + res_off)
        np.savetxt('contrast-cw-new.txt', np.vstack([scheduler._freqs, contrasts]))
        # # plot_freq_contrast(freqs, contrasts, fname='contrast-by-single')


if __name__ == '__main__':
    # unittest.main()

    st = SingleFrequencyTest()
    # st.test_single_freq(1)

    # st.test_single_freq(10, 2.87)
    # print('==='*10)
    # st.test_single_freq(10, 2.85)

    st.test_scanning()
