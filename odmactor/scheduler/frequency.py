"""
ODMR detecting classes
---
1. Continuous-wave detection (frequency-domain method)
Basic: 激光和微波同时施加，单光子探测器持续探测全过 程中的光子数，通过扫描微波频 率产生频谱。
Remarks: 连续波谱不是典型 的量子传感过程，实验中系统处 于开放稳态，相干性在此实验中几乎不起作用
---
2. Pulse detection (frequency-domain method)
输入：
	- 微波频率范围和步长（起止点/中心和宽度）--> freqs = linspace(,,)
	- ASG 序列周期个数，serial time intervals parameters. e.g. 10ms * 100 = 1s ，后者作为微波在单个频率点的持续时间 T = N * t
	- 单个周期内激光初始化时间 t_init, 微波翻转时间 t_mw，读出时间 t_read
	-
输出：对比度数据，
调度过程  结果读出 ：
	- 生成两通道数据（channel['laser']，channel['mw'] --> asg_data: Matrix）
	- For each freq in freqs: Asg start --> Asg stop 历时 T ，存储 APD 数据 到 self._cache
	- 最后统一计算各频率点的对比度，结果到 self.result : {'freqs': …, 'ratio': …}，self.result_detail
经验参数：
    - init: 5us
    - init-mw-interval 3us
    - mw: 5us
    - readout sign: 400ns
    - readout-interval: 200ns
    - readout ref: 400s
"""

import threading
import time
import numpy as np
import scipy.constants as C
from typing import List
from odmactor.scheduler.base import FrequencyDomainScheduler
from odmactor import utils


class CWScheduler(FrequencyDomainScheduler):
    """
    Continuous-Wave ODMR scheduler
    """

    def __init__(self, *args, **kwargs):
        super(CWScheduler, self).__init__(*args, **kwargs)
        self.name = 'CW ODMR Scheduler'

    def configure_odmr_seq(self, t, N: int):
        """
        Wave form for single period:
            laser (no asg control sequence):
            ---------------------------------
            |                               |
            |                               |
            microwave (no asg control sequence):
            ---------------------------------
            |                               |
            |                               |
            asg tagger acquisition channel (continuous readout in fact):
            ---------------------------------
            |                               |
            |                               |
        All units for the parameters is 'ns'
        :param t: binwidth parameter for TimeTagger.Counter
        :param N: n_values for TimeTagger.Counter
        """
        cont_seq = [t, 0]
        if self.mw_ttl == 0:
            mw_seq = utils.flip_sequence(cont_seq)
        else:
            mw_seq = cont_seq
        self.download_asg_sequences(laser_seq=cont_seq, mw_seq=mw_seq, tagger_seq=cont_seq, N=N)

    def run_single_step(self, power, freq, mw_control='on') -> List[float]:
        """
        Single-frequency & single-power setting for running the scheduler
        :param power: MW power, unit: dBm
        :param freq: MW frequency, unit: Hz
        :param mw_control: 'on' or 'off'
        :return: 1-D array: [N,]
        """
        print('running for freq = {:.3f} GHz, power = {:.2f} dBm ...'.format(freq / C.giga, power))

        # set MW parameters
        self.configure_mw_paras(power, freq)

        mw_seq_on = self.mw_control_seq()
        if mw_control == 'off':
            self.mw_control_seq([0, 0])
        elif mw_control == 'on':
            pass
        else:
            raise ValueError('unsupported mw_control parameter')

        self._start_device()
        time.sleep(0.5)  # 先让激光和微波开几秒
        time.sleep(self.asg_dwell)
        counts = self.counter.getData().ravel().tolist()
        self.stop()

        # recover the asg control sequence for MW to be 'on'
        if mw_control == 'off':
            self.mw_control_seq(mw_seq_on)

        return counts

    def _scan_freqs_and_get_data(self, with_ref: bool = True):
        """
        Override this function in CWScheduler class
        :param with_ref: if True, detect both signals and reference signals in two sequent asg_dwell periods
        """
        if with_ref:
            for i, freq in enumerate(self._freqs):
                self._mw_instr.write_float('FREQUENCY', freq)
                print('scanning freq {:.3f} GHz'.format(freq / C.giga))
                #################################
                # 1. MW on
                self._mw_instr.write_bool('OUTPUT:STATE', True)
                # print('MW on/off status:', self._mw_instr.instrument_status_checking)

                t = threading.Thread(target=self._get_data, name='thread-on-{}'.format(i))
                time.sleep(self.time_pad)
                time.sleep(self.asg_dwell)  # accumulate counts
                t.start()  # begin readout
                time.sleep(self.time_pad)

                # 2. MW off
                self._mw_instr.write_bool('OUTPUT:STATE', False)
                # print('MW on/off status:', self._mw_instr.instrument_status_checking)
                t = threading.Thread(target=self._get_data_ref, name='thread-off-{}'.format(i))
                time.sleep(self.time_pad)
                time.sleep(self.asg_dwell)  # accumulate counts
                t.start()  # begin readout
                time.sleep(self.time_pad)

            print('finished data acquisition')
        else:
            super(CWScheduler, self)._scan_freqs_and_get_data()

    def _cal_counts_result(self, with_ref: bool = False):
        """
        Override this function in CWScheduler class
        :param with_ref: if True, detect both signals and reference signals in two sequent asg_dwell periods
        """
        super(CWScheduler, self)._cal_counts_result()
        if with_ref:
            counts_ref = [np.mean(ls) for ls in self._data_ref]
            self._result.append(counts_ref)
            self._result_detail['counts_ref'] = counts_ref
            self._result_detail['origin_data_ref'] = self._data_ref

    def _acquire_data(self, with_ref: bool = False):
        """
        Scanning time intervals to acquire data for Time-domain Scheduler
        """
        # 1. scan time intervals
        self._scan_freqs_and_get_data(with_ref)

        # 2. calculate result (count with/without reference)
        self._cal_counts_result(with_ref)

        # 3. save result
        self.save_result(self._gene_data_result_fname())

    def run_scanning(self, mw_control='on', with_ref: bool = False):
        """
        Run the scheduler under scanning-frequency mode
        1) start device
        2) acquire data timely
        :param with_ref: if True, detect both signals and reference signals in two sequent asg_dwell periods
        """
        mw_seq_on = self.mw_control_seq()
        if mw_control == 'off':
            self.mw_control_seq([0, 0])
        elif mw_control == 'on':
            pass
        else:
            raise ValueError('unsupported MW control parameter (should be "on" or "off"')

        print('Begin to run {}. Frequency: {:.3f} - {:.3f} GHz.'.format(self.name, self._freqs[0] / C.giga,
                                                                        self._freqs[-1] / C.giga))
        print('t: {:.2f} ns, N: {}, T: {:.2f} s, n_freqs: {}'.format(self._asg_conf['t'] / C.nano, self._asg_conf['N'],
                                                                     self.mw_dwell, len(self._freqs)))
        print('Estimated total running time: {:.2f} s'.format(self.time_total))

        self._data_ref.clear()  # this list is used only in CWScheduler
        self._start_device()
        self._acquire_data(with_ref)  # scanning MW frequencies in this loop
        self.stop()

        # 恢复微波的ASG的MW通道为 on
        if mw_control == 'off':
            self._asg_sequences[self.channel['mw'] - 1] = mw_seq_on
            self.asg_connect_and_download_data(self._asg_sequences)

    # def run_origin(self, mw_control='on'):
    #     """
    #     1) start device
    #     2) acquire data timely
    #     """
    #     mw_seq_on = self._asg_sequences[self.channel['mw'] - 1]
    #     if mw_control == 'off':
    #         mw_seq_off = [0, sum(mw_seq_on)]
    #         self._asg_sequences[self.channel['mw'] - 1] = mw_seq_off
    #         self.asg_connect_and_download_data(self._asg_sequences)
    #     elif mw_control == 'on':
    #         pass
    #     else:
    #         raise ValueError('unsupported MW control parameter (should be "or" or "off"')
    #
    #     self._start_device()
    #     self._acquire_data()
    #
    #     # 恢复微波的ASG的MW通道为 on
    #     self._asg_sequences[self.channel['mw'] - 1] = mw_seq_on
    #     self.asg_connect_and_download_data(self._asg_sequences)
    #
    #     # self._start_device()
    #     # self._acquire_data()
    #
    #     data_sig = []
    #     data_ref = []
    #     t = self._asg_conf['t']  # unit: s
    #     num_freqs = len(self._freqs)
    #     self.counter = tt.Counter(self.tagger, channels=[1], binwidth=int(self._asg_conf['t'] / C.pico),
    #                               n_values=self._asg_conf['N'])
    #
    #     for freq in self._freqs:
    #         print('scanning freq {:.4f} GHz'.format(freq / C.giga))
    #         # self._mw_instr.write_bool('OUTPUT:STATE', True)
    #         self._mw_instr.write_float('FREQUENCY', freq)
    #         time.sleep(0.2)
    #
    #         # self.laser_on_seq(t / C.nano)
    #         self.mw_on_seq(t / C.nano)
    #         time.sleep(self.asg_dwell)
    #         data_sig.append(self.counter.getData())
    #         self.counter.clear()
    #         # print('---------- on')
    #         # print(self.counter.getData().ravel())
    #         # print(self._asg_sequences)
    #         # self.counter.clear()
    #         # print('after clear:', self.counter.getData())
    #
    #         self.asg.stop()  # stop firstly
    #         # self.laser_off_seq()
    #         self.mw_off_seq()
    #         time.sleep(self.asg_dwell)
    #         data_ref.append(self.counter.getData())
    #         self.counter.clear()
    #         # print('---------- off')
    #         # print(self.counter.getData().ravel())
    #         # print(self._asg_sequences)
    #         # self.counter.clear()
    #         # print('after clear:', self.counter.getData())
    #
    #     self.asg.stop()
    #     self.counter.stop()
    #
    #     data_sig_avg = [np.mean(ls) for ls in data_sig]
    #     data_ref_avg = [np.mean(ls) for ls in data_ref]
    #     contrast = [np.abs(data_sig_avg[i] - data_ref_avg[i]) / data_ref_avg[i] for i in range(num_freqs)]
    #
    #     self._result = [self._freqs, contrast]
    #
    #     self._result_detail = {
    #         'freqs': self._freqs,
    #         'contrast': contrast,
    #         'data_ref': data_sig,
    #         'data_sig': data_ref
    #     }
    #
    #     fname = os.path.join(self.output_dir,
    #                          'CW-ODMR-result-{}-{}'.format(str(datetime.date.today()), round(time.time())))
    #     with open(fname + '.json', 'w') as f:
    #         json.dump(self._result_detail, f)
    #     print('data has been saved into {}'.format(fname))


class PulseScheduler(FrequencyDomainScheduler):
    """
    Pulse-based ODMR manipulation scheduler
    """

    def __init__(self, *args, **kwargs):
        super(PulseScheduler, self).__init__(*args, **kwargs)
        self.name = 'Pulse ODMR Scheduler'

    def set_odmr_seq(self, laser_seq: list = None, mw_seq: list = None, tagger_seq: list = None):
        if self.mw_ttl == 0:
            mw_seq = utils.flip_sequence(mw_seq)
        self.download_asg_sequences()

    def configure_odmr_seq(self, t_init, t_mw, t_read_sig, t_read_ref=None, inter_init_mw=3000, pre_read=200,
                           inter_mw_read=500, inter_readout=200, inter_period=200, N: int = 1000):
        """
        Wave form for single period:
            asg laser channel:
            -----               -------------
            |   |               |           |
            |   |---------------|           |
            asg microwave channel:
                    -------------
                    |           |
            --------|           |------------
            asg tagger acquisition channel:
                                ---   ---
                                | |   | |
            --------------------| |---| |----
        All units for the parameters is 'ns'
        :param t_init: time span for laser initialization, e.g. 5000
        :param t_mw: time span for microwave actual operation in a ASG period, e.g. 800
        :param t_read_sig: time span for fluorescence signal readout, e.g. 400
        :param t_read_ref: time span for reference signal readout, e.g. 400
                            if the parameter is not assigned, means single-pulse readout
        :param inter_init_mw: time interval between laser initialization and MW operation pulses, e.g. 3000
        :param inter_readout: interval between single readout pulse and reference signal readout, e.g. 200
                            when t_read_ref is assigned, this parameter will play its role
        :param inter_period: interval between two neighbor periods, e.g. 200
        :param N: number of ASG operation periods
        """
        # unit: ns
        # total time for 'N' period, also for MW operation time at each frequency point
        if t_read_ref is not None:
            self.two_pulse_readout = True
            laser_seq = [t_init, inter_init_mw + t_mw + inter_mw_read,
                         pre_read + t_read_sig + inter_readout + t_read_ref + inter_period, 0]
            mw_seq = [0, t_init + inter_init_mw, t_mw,
                      inter_mw_read + pre_read + t_read_sig + inter_readout + t_read_ref + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw + inter_mw_read + pre_read, t_read_sig, inter_readout,
                          t_read_ref, inter_period]
            # apd_seq = [sum(tagger_seq[:-4]), sum(tagger_seq[-4:])]
        else:
            # single-pulse readout
            laser_seq = [t_init, inter_init_mw + t_mw + inter_mw_read, pre_read + t_read_sig + inter_period, 0]
            mw_seq = [0, t_init + inter_init_mw, t_mw, inter_mw_read + pre_read + t_read_sig + inter_period]
            tagger_seq = [0, t_init + inter_init_mw + t_mw + inter_mw_read + pre_read, t_read_sig, inter_period]
            # apd_seq = [sum(tagger_seq[:-2], sum(tagger_seq[-2:]))]

        if self.mw_ttl == 0:
            #     # low-level effective
            mw_seq = utils.flip_sequence(mw_seq)

        self.download_asg_sequences(laser_seq, mw_seq, tagger_seq, N)

    def run_single_step(self, freq, power=None, mw_control='on') -> List[float]:
        """
        Single-frequency & single-power setting for running the scheduler
        :param freq: MW frequency, unit: Hz
        :param power: MW power, unit: dBm
        :param mw_control: 'on' or 'off'
        :return: 1-D array: [N,]
        """
        print('running for freq = {:.3f} GHz, power = {:.2f} dBm ...'.format(freq / C.giga, power))

        # set MW parameters
        self.configure_mw_paras(power, freq)

        # start sequence for time: N*t
        mw_seq_on = self.mw_control_seq()
        if mw_control == 'off':
            self.mw_control_seq([0, 0])
        elif mw_control == 'on':
            pass
        else:
            raise ValueError('unsupported mw_control parameter')

        self._start_device()
        time.sleep(0.5)  # 先让激光和微波开几秒
        time.sleep(self.asg_dwell)
        counts = self.counter.getData().ravel().tolist()
        self.stop()

        # recover the asg control sequence for MW to be 'on'
        if mw_control == 'off':
            self.mw_control_seq(mw_seq_on)

        return counts
