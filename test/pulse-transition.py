from odmactor.scheduler import PulseScheduler
import scipy.constants as C
import time
import matplotlib.pyplot as plt
import unittest
from odmactor.utils.plotting import plot_freq_contrast

"""
Steps:
    - MW
    - ASG sequences
    - counter
    - run
"""

freq_start = 2.85 * C.giga
freq_end = 2.89 * C.giga
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
tagger_input = {'apd': 1, 'asg': 2}

t_init = 2e4
t_mw = 5e4
inter_init_mw = 1e4

transition_time = 0.0  # 平衡过渡时间


def contrast_testing():
    scheduler = PulseScheduler(transition_time=transition_time)
    scheduler.channel = channel_dict
    scheduler.tagger_input = tagger_input
    scheduler.configure_mw_paras(power=10)
    scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)
    scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)
    scheduler.configure_tagger_counting(reader='cbm')
    scheduler.run()
    scheduler.close()
    with open('seq.txt', 'w') as f:
        f.write(scheduler.sequences_strings)
    # scheduler.sequences_figure.save('pulse-seq.png', dpi=400)

    plot_freq_contrast(*scheduler.result, fname='contrasts-two-pulse')

def counts_testing():
    scheduler = PulseScheduler(transition_time=transition_time)
    scheduler.channel = channel_dict
    scheduler.tagger_input = tagger_input
    scheduler.configure_mw_paras(power=10)
    scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)
    scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)
    scheduler.configure_tagger_counting()

    scheduler.run('on')
    scheduler.close()
    res_on = scheduler.result

    # scheduler.run('off')
    # scheduler.close()
    # res_off = scheduler.result
    #
    # contrasts = [abs(c_off - c_on) / c_on for c_on, c_off in zip(res_on[1], res_off[1])]

    plt.plot(*res_on, label='MW on')
    # plt.plot(*res_off, label='MW off')
    plt.title('Counts comparison')
    plt.savefig('counts_on_off', dpi=400)

    # plot_freq_contrast(res_on[0], contrasts, fname='contrasts-two-pulse')

#
# class PulseTest(unittest.TestCase):
#     # def test_scanning_freqs_single_readout(self):
#     #     scheduler = PulseScheduler(transition_time=transition_time)
#     #     scheduler.channel = channel_dict
#     #     scheduler.tagger_input = tagger_input
#     #     scheduler.configure_mw_paras(power=10)
#     #     scheduler.configure_odmr_seq(t_init, t_mw, inter_init_mw, N=N)
#     #     print('当前序列：')
#     #     print(scheduler.sequences_strings)
#     #     print()
#     #     scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)
#     #     scheduler.configure_tagger_counting()
#     #
#     #     scheduler.run('on')
#     #     scheduler.close()
#     #     res_on = scheduler.result
#     #
#     #     scheduler.run('off')
#     #     scheduler.close()
#     #     res_off = scheduler.result
#     #
#     #     contrasts = [abs(c_off - c_on) / c_on for c_on, c_off in zip(res_on[1], res_off[1])]
#     #
#     #     plt.plot(*res_on, label='MW on')
#     #     plt.plot(*res_off, label='MW off')
#     #     plt.title('Counts comparison')
#     #     plt.savefig('counts_on_off', dpi=400)
#     #
#     #     plot_freq_contrast(res_on[0], contrasts, fname='contrasts-two-pulse')
#
#     def test_scanning_freqs_two_pulse_readout(self):
#         scheduler = PulseScheduler(transition_time=transition_time)
#         scheduler.channel = channel_dict
#         scheduler.tagger_input = tagger_input
#         scheduler.configure_mw_paras(power=10)
#         scheduler.configure_odmr_seq(t_init, t_mw, inter_init_mw, N=N)
#         scheduler.configure_tagger_counting(reader='cbm')
#         scheduler.run()
#         scheduler.close()
#
#         plot_freq_contrast(*scheduler.result, fname='contrasts-two-pulse')
#

# if __name__ == '__main__':
#     unittest.main()
def wave_form():
    scheduler = PulseScheduler(transition_time=transition_time)
    scheduler.channel = channel_dict
    scheduler.tagger_input = tagger_input
    scheduler.configure_mw_paras(power=10)
    scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)
    scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)
    scheduler.configure_tagger_counting()

    scheduler.sequences_figure.savefig('seq.png', dpi=400)

if __name__ == '__main__':
    counts_testing()