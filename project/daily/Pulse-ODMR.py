from scheduler.pulse import PulseScheduler
import scipy.constants as C
import time
from utils.plotting import plot_freq_contrast

freq_start = 2.8 * C.giga
freq_end = 2.9 * C.giga
freq_step = 2 * C.mega

##########################
# e.g.
# N = 100 000
# single ASG period: 5 us
# MW period at each freq: 5*N ~ 0.5s
# total time: 5 * N * 10(11) ~ 5 s
##########################
N = int(1e4)
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


scheduler = PulseScheduler()

scheduler.channel = channel_dict
scheduler.tagger_input = tagger_input
scheduler.configure_odmr_seq(t_mw=200000, t_init=2000, t_read_sig=800, t_read_ref=800, t_interval=50,
                             N=N)  # unit: ns
# scheduler.configure_odmr_seq(t_mw=200000000, t_init=100000000, t_read_sig=100000000, t_read_ref=100000000,
#                              t_interval=100000000,
#                              N=N)  # unit: ns
scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)
scheduler.configure_mw_paras(power=25)

print('单个频率周期：{:.2f} s，总时间:{:.2f} s'.format(scheduler.mw_dwell, scheduler.time_total))


scheduler.run()  # auto stop
scheduler.close()

print('sync delay: {} ns'.format(scheduler.sync_delay))
print('total running time required: {:.4f} s'.format(scheduler.time_total))
# time.sleep(scheduler.time_total + 10)  # main thread sleeps for time_total


# print(scheduler.result.round(2))
print(scheduler.result)
print()
print(scheduler.result_detail)

res = scheduler.result

print()
print(res[0])
print(res[1])
# print(res)
plot_freq_contrast(scheduler.result[0], scheduler.result[1])
