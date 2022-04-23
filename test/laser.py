import scipy.constants as C
import time
from odmactor.scheduler import PulseScheduler

"""
Pulse-ODMR 的序列来测试激光是否被正确控制
"""

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
tagger_input = {
    'apd': 1,
    'asg': 2
}

scheduler = PulseScheduler()

scheduler.channel = channel_dict
scheduler.tagger_input = tagger_input

t_mw = 2e4
t_init = 1e4
t_read_sig = 1e4

t_read_ref = 1e4
t_interval = 5e3

scheduler.configure_odmr_seq(t_mw=t_mw, t_init=t_init, t_read_sig=t_read_sig, t_read_ref=t_read_ref,
                             t_interval=t_interval, N=N)  # unit: ns
# scheduler.asg_data = [
#     [t_init, t_mw + 2 * t_interval, t_read_sig + t_interval + t_read_ref, t_interval],
#     [0, t_init + t_interval, t_mw, 3 * t_interval + t_read_sig + t_read_ref],
#     [0, 0],
#     [0, 0],
#     [0, t_mw + 2 * t_interval + t_init, t_read_sig, t_interval, t_read_ref, t_interval],
#     [0, 0],
#     [0, 0],
#     [0, 0],
# ]


scheduler.asg.start()

time.sleep(30)
scheduler.asg.stop()

scheduler.asg.close_device()
