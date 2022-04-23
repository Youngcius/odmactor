import numpy as np

from odmactor.scheduler import CWScheduler
import scipy.constants as C


def scan_freqs_CW(start, end, step, t, N: int):
    scheduler = CWScheduler()
    scheduler.configure_mw_paras(power=10)
    scheduler.configure_odmr_seq(t, N)
    scheduler.set_mw_scan_freq_start_stop(start, end, step)
    scheduler.configure_tagger_counting()  # MW always on or off
    scheduler.run()

    res = scheduler.result  # [freqs, counts]
    np.savetxt('scan-cw-f-c.txt', res)



if __name__ == '__main__':
    t_ns = 1e5
    N = int(1e4)
    freq_start = 2.85 * C.giga
    freq_end = 2.89 * C.giga
    freq_step = 1 * C.mega

    scan_freqs_CW(freq_start, freq_end, freq_step, t_ns, N)
