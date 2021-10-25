from scheduler.continuity import CWScheduler
import scipy.constants as C
from utils.plotting import plot_freq_contrast

scheduler = CWScheduler()

binwidth = 400000

N = int(1e4)  # 每个频率重复测量 N 次
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

freq_start = 2.85 * C.giga
freq_end = 2.9 * C.giga
freq_step = 5 * C.mega

scheduler.channel = channel_dict
scheduler.tagger_input = tagger_input
scheduler.configure_mw_paras(power=20)  # 20dBm, 3GHz
scheduler.configure_odmr_seq(t=binwidth, N=N)
# scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)  # dwell=None
scheduler.set_mw_freqs(freq_start, freq_end, freq_step)
# print(scheduler._asg_conf)


print('单个频率周期：{:.2f} s，总时间:{:.2f} s'.format(scheduler.mw_dwell, scheduler.time_total))



# print(scheduler._asg_sequences)


scheduler.run()  # auto stop
scheduler.close()

# print(scheduler.result.round(2))
print(scheduler.result)
print(scheduler.result_detail)

res = scheduler.result

print(res)
plot_freq_contrast(res[0], res[1])

