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

scheduler = PulseScheduler()

scheduler.channel = channel_dict
scheduler.configure_odmr_seq(t_mw=2000, t_init=100, t_read_sig=500, t_read_ref=500, t_interval=50,
                             N=N)  # unit: ns
# scheduler.configure_odmr_seq(t_mw=200000000, t_init=100000000, t_read_sig=100000000, t_read_ref=100000000,
#                              t_interval=100000000,
#                              N=N)  # unit: ns
scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)  # dwell=None
scheduler.configure_tagger_counting()
scheduler.run()

print('sync delay: {} ns'.format(scheduler.sync_delay))
print('total running time required: {:.4f} s'.format(scheduler.time_total))
time.sleep(scheduler.time_total + 10)  # main thread sleeps for time_total
# scheduler.stop()
scheduler.close()
