import matplotlib.pyplot as plt
from odmactor.scheduler import CWScheduler, PulseScheduler
import scipy.constants as C
import numpy as np
import time
from odmactor.instrument import *


# 50 Hz

f = 50
period = 1 / f / C.nano

# t_ns = 1e5
N = int(1e3 / 4)


freq_start = 0.04 * C.giga
freq_end = 0.07 * C.giga
freq_step = 10 * C.mega
p = -15

# %%
scheduler = CWScheduler(use_lockin=True, mw_ttl=1, with_ref=True)
scheduler.reconnect()
scheduler.configure_mw_paras(p)
scheduler.configure_odmr_seq(period, N)
scheduler.configure_lockin_counting()
scheduler.set_mw_freqs(freq_start, freq_end, freq_step)
# scheduler.configure_tagger_counting()

scheduler.run_scanning()
scheduler.close()


print(scheduler._data)

print(scheduler._data_ref)


