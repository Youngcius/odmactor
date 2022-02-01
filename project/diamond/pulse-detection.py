"""
Create date: 2022-01-01
Test object: diamond bulk
"""

from odmactor.scheduler import PulseScheduler
import scipy.constants as C
import time
import matplotlib.pyplot as plt
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
freq_step = 1 * C.mega

##########################
# e.g.
# N = 100 000
# single ASG period: 5 us
# MW period at each freq: 5*N ~ 0.5s
# total time: 5 * N * 10(11) ~ 5 s
##########################
N = int(1e5)
channel_dict = {
    'laser': 1,
    'mw': 2,
    'apd': 3,
    'tagger': 5
}
tagger_input = {'apd': 1, 'asg': 2}

t_init = 3e3
t_mw = 800
inter_init_mw = 5e3


# transition_time = 0.0  # 平衡过渡时间


def contrast_testing():
    scheduler = PulseScheduler()
    scheduler.channel = channel_dict
    scheduler.tagger_input = tagger_input
    scheduler.configure_mw_paras(power=10)
    scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)
    scheduler.set_mw_freqs(freq_start, freq_end, freq_step)
    scheduler.configure_tagger_counting(reader='cbm')
    scheduler.run()
    scheduler.save_result('pulse-result-1123')
    scheduler.close()
    # with open('seq.txt', 'w') as f:
    #     f.write(scheduler.sequences_strings)
    scheduler.sequences_figure.savefig('pulse-seq.png', dpi=400)

    plot_freq_contrast(*scheduler.result, fname='contrasts-two-pulse')


def single_freq_testing(freq):
    scheduler = PulseScheduler()
    scheduler.channel = channel_dict
    scheduler.tagger_input = tagger_input
    scheduler.configure_mw_paras(power=10, freq=freq)
    scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)
    scheduler.configure_tagger_counting(reader='cbm')
    scheduler.run(scan=False)
    scheduler.save_result('pulse-result-single-freq')
    scheduler.close()
    # with open('seq.txt', 'w') as f:
    #     f.write(scheduler.sequences_strings)

    # plot_freq_contrast(*scheduler.result, fname='contrasts-two-pulse')


def counts_testing():
    scheduler = PulseScheduler()
    scheduler.channel = channel_dict
    scheduler.tagger_input = tagger_input
    scheduler.configure_mw_paras(power=10)
    scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, inter_init_mw=inter_init_mw, N=N)
    scheduler.sequences_figure.savefig('pulse-seq-test.png', dpi=400)

    scheduler.set_mw_freqs(freq_start, freq_end, freq_step)
    scheduler.configure_tagger_counting(reader='cbm')

    scheduler.run_scanning('on')
    # run之后自动stop各个仪器，
    # scheduler.close()

    res_on = scheduler.result

    scheduler.run_scanning('off')
    # 自动重启各个仪器
    scheduler.close()
    res_off = scheduler.result

    contrasts = [abs(c_off - c_on) / c_on for c_on, c_off in zip(res_on[1], res_off[1])]

    plt.plot(*res_on, label='MW on')
    plt.plot(*res_off, label='MW off')
    plt.title('Counts comparison')
    plt.savefig('counts_on_off (Pulse)', dpi=400)

    plot_freq_contrast(res_on[0], contrasts, fname='contrasts-single-pulse-repetition')


def wave_form():
    scheduler = PulseScheduler()
    scheduler.channel = channel_dict
    scheduler.tagger_input = tagger_input
    scheduler.configure_mw_paras(power=10)
    scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)
    scheduler.set_mw_freqs(freq_start, freq_end, freq_step)
    scheduler.configure_tagger_counting()

    scheduler.sequences_figure.savefig('seq.png', dpi=400)
    scheduler.close()


if __name__ == '__main__':
    contrast_testing()
    # wave_form()
    # counts_testing()
    # single_freq_testing(2.85e9)
    # data = np.loadtxt('data.txt')
    # sig = data[::2]
    # ref = data[1::2]
    #
    # plt.hist(sig, bins=15)
    # plt.hist(ref, bins=15)
    # plt.show()
    # sigavg = np.mean(sig)
    # refavg = np.mean(ref)
    # print(sigavg, refavg)
