from scheduler.continuity import CWScheduler
import scipy.constants as C
from utils.plotting import plot_freq_contrast

scheduler = CWScheduler()

N = int(4e6) # 每个频率重复测量 N 次
channel_dict = {
    'laser': 1,
    'mw': 2,
    'apd': 3,
    'tagger': 5
}
tagger_input = {
    'apd': 1,
    'asg': 2
}

freq_start = 2 * C.giga
freq_end = 2.87 * C.giga
freq_step = 0.87 * C.mega

# scheduler.time_total = 1
scheduler.channel = channel_dict
scheduler.tagger_input = tagger_input
scheduler.configure_mw_paras(power=10)  # 20dBm, 3GHz
scheduler.configure_odmr_seq(t=250, t_read=200, N=N)
scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)  # dwell=None

scheduler.run()  # auto stop
print('单个频率周期：{}，总时间:{}'.format(scheduler.mw_dwell, scheduler.time_total))
scheduler.close()

print(scheduler.result.round(2))
plot_freq_contrast(scheduler.result[0], scheduler.result[1])
