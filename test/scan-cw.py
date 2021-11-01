from odmactor.scheduler import CWScheduler
import scipy.constants as C

freq_start = 2.85 * C.giga
freq_end = 2.89 * C.giga
freq_step = 1 * C.mega


def scan_freqs_CW(start, end, step, t, N:int):
    scheduler = CWScheduler()
    scheduler.configure_mw_paras(power=10)
    scheduler.configure_tagger_counting() # MW always on or off
    scheduler.run()
