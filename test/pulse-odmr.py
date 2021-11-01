from odmactor.scheduler.pulse import PulseScheduler
import scipy.constants as C
import time

freq_start = 2.5 * C.giga
freq_end = 3.5 * C.giga
freq_step = 100 * C.mega

##########################
# e.g.
# N = 100 000
# single ASG period: 5 us
# MW period at each freq: 5*N ~ 0.5s
# total time: 5 * N * 10(11) ~ 5 s
##########################
N = int(1e6)
channel_dict = {
    'laser': 1,
    'mw': 2,
    'apd': 3,
    'tagger': 5
}

t_init = 2e4
t_mw = 5e4
inter_init_mw = 1e4


scheduler = PulseScheduler()

scheduler.channel = channel_dict
scheduler.configure_mw_paras(power=10)
scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)

scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)  # dwell=None
scheduler.configure_tagger_counting()
scheduler.run()

# print('sync delay: {} ns'.format(scheduler.sync_delay))
# print('total running time required: {:.4f} s'.format(scheduler.time_total))
# time.sleep(scheduler.time_total + 10)  # main thread sleeps for time_total
# scheduler.stop()
scheduler.close()


