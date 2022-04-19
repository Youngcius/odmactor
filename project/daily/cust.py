from odmactor.scheduler import FrequencyDomainScheduler


class MyScheduler(FrequencyDomainScheduler):
    def __init__(self, *args, **kwargs):
        super(MyScheduler, self).__init__(*args, **kwargs)
        self.name = 'My Scheduler'

    def configure_odmr_seq(self, t_read_sig, t_read_ref, pre_read, inter_read, inter_period, N=100000):
        tagger_seq = [0, pre_read, t_read_sig, inter_read, t_read_ref, inter_period]
        t = sum(tagger_seq)
        laser_seq = [t, 0]

        # configure ASG period information
        self._conf_time_paras(t, N)

        idx_laser_channel = self.channel['laser'] - 1
        idx_tagger_channel = self.channel['tagger'] - 1

        self.reset_asg_sequence()
        self._asg_sequences[idx_laser_channel] = laser_seq
        self._asg_sequences[idx_tagger_channel] = tagger_seq

        # connect & download pulse data
        self.asg_connect_and_download_data(self._asg_sequences)

    def run_scanning(self, mw_control: str = 'on'):
        pass


if __name__ == '__main__':

    scheduler = MyScheduler()
    # scheduler.set_mw_freqs()
    scheduler.run_scanning()