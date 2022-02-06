# Odmactor (ODMR Actor)

ODMR software SDK integrating functions of ODMR detection and spin manipulation.

*This project is on continuous updating & using ...*

> Last updation date: February, 2022<br>
> Current implementation: ODMR detection

## Design Method

### Measurement modes and scheduling strategies

1. Frequency-domain detection<br>
   *Scan MW frequencies*
    - Continuous-Wave ODMR
        1) polarize spin systems
        2) operate continuous MW and readout signals for the whole sequence period
    - Pulse ODMR
        1) polarize spin systems
        2) apply a fixed-interval MW pulse
        3) readout final fluorescence signals
2. Time-domain detection<br>
   *Scan time intervals*
    - Ramsey detection
        1) initialize spin to ground state
        2) initialize spin to equal-amplitude superposition state using calibrated MW $\frac{\pi}{2}$ pulse
        3) wait for a time interval $\tau$
        4) operate a calibrated MW $\frac{\pi}{2}$ pulse again
        5) readout final spin state (population)
    - Rabi oscillation
        1) initialize spin to ground state
        2) operate a short MW pulse (100~200 ns)
        3) readout final spin state (population)
    - T1 relaxation
        1) initialize spin to excited state using calibrated MW $\pi$ pulse
        2) wait a time interval $\tau$
        3) readout final spin state (population)


3. Dynamical Decoupling<br>
   *Currently not implemented*
    - Hahn echo
    - High-order Dynamical Decoupling
4. Spin manipulation<br>
   *Currently not implemented*

The above specific scheduling methods could be abstracted into different detection sequences in experiments. They are
all controlled in precision of "ns".

![odmr-sequence](./asset/figure/odmr-sequence.png)

Thus the series of ODMR measurement experimetns is simplified by a "pipeline". CW and Pulse ODMR could be used to
characterize environmental physical quantities, while they also could be used to calibrate energy gap or fine
structures. Ramsey detection is usually used to characterize T2* (dephasing time) of spin systems. Some typical results
of them are like the following figure.

![](./asset/figure/odmr-magnet.png)

### OOP Classes implementation

![Scheduler classes profile](./asset/figure/classes.png)

**general scheduling methods**

- `scheduler.config_odmr_seq()`: configure ASG control sequences for laser, MW and tagger
- `schedulel.stop()`: stop all hardware (ASG, MW, Tagger) scheduling
- `schedulel.close()`: release instrument (ASG, MW, Tagger) resources
- `schedulel.save_result()`: save detailed measurement result into a ".json" file

**necessary data fields of schedulers**

- `scheduler.channel`: a `dict` instance in Python, presenting ASG channels to output control sequences
- `scheduler.tagger_input`: a `dict` instance in Python, presenting Tagger input channels (signal, trigger, etc.)
- `scheduler.result`: a `list` instance in Python, consisting of aggregated (avg or mean) measurement result
- `scheduler.result_detail`: a `dict` instance in Python, consisting of detailed measurement result
- `scheduler.pi_pulse`: `{'freq': ..., 'power': ..., 'time': ...}`, consisting of configuration information of the
  calibrated MW $\pi$ pulse
- `scheduler.sequences_figure`: visualized ASG control sequences

**specific scheduling methods**

- `scheduler.run_scanning()`: for `FrequencyDomainScheduler` and `TimeDomainScheduler`, scan time intervals and MW
  frequencies, respectively
- `scheduler.run_single_step()`: for `FrequencyDomainScheduler` and `TimeDomainScheduler`,
- `scheduler.set_delay_times()`: for `TimeDomainScheduler`, this function should be called to design time interval
  scanning points
- `scheduler.set_mw_freqs()`: for `FrequencyDomainScheduler`, this function should be called to design MW frequency
  scanning points

## Hardware support

### General ODMR platform

![ODMR platform](./asset/figure/platform.png)

### Specific hardware used

Most general hardware resources are supported by our Odmactor programs, while some specific instrument currently used in
our [QIM](https://quantum.lab.arizona.edu/) group are as follows.

**Arbitrary Sequence Generator (ASG)**

The most important instrument to synchronize control measurement processes.

- Vendor: [CIQTEK](https://www.ciqtek.com/)
- Type: ASG8005
- Output: 8 control channels

**Microwave instrument (MW)**

- Vendor:
- Type:
- Frequency range:

**Time Tagger (Tagger)**

This is a T/D convertor as well as A/D convertor necessary for data acquisition.

- Vendor: [Swabian](https://www.swabianinstruments.com/)
- Type: Time Tagger 20

## Implementation method

### Basic Implementation Steps

- Configuration
    1. Set channels of ASG
    2. Set parameters (e.g. frequencies, power) of MW
    3. Configure pulse sequences for ODMR experiments
    4. Configure the counting module of Tagger
- Start devices & Acquire data
    1. run the "scheduler" (e.g. execute `scheduler.run()` method)
    2. the result will be saved into properties `scheduler.result` and `scheduer.result_detail`
- Analysis & Visualization
    1. visualize your result (e.g. contrast figure, counting figure corresponding to frequencies)
    2. use suitable analytical functions to fit your data curve

### Code example

More examples could be found in the "./project/tutorial/" folder. Herein is a typical scheduling pipeline instance.

```python
# ------- Pulse ODMR Measurement -------
# configure parameters
channel_dict = {
    'laser': 1,
    'mw': 2,
    'apd': 3,
    'tagger': 5
}
tagger_input = {'apd': 1, 'asg': 2}

scheduler = PulseScheduler()  # Pulse-ODMR scheduler instance
scheduler.channel = channel_dict  # set ASG control channels
scheduler.tagger_input = tagger_input  # set tagger input channels
scheduler.configure_mw_paras(power=10) # configure MW power
scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)
scheduler.set_mw_freqs(freq_start, freq_end, freq_step) # set measured frequencies
scheduler.configure_tagger_counting(reader='cbm')  # 'Counter Between Markers' measurement mode

# run
scheduler.run_scanning()
scheduler.save_result(fname) # actually, this method is automatically called in the run_scanning() method

# result analysis
counts_sig_ref = scheduler.result # [freqs, counts, counts_reference]
contrast = [sig / ref for sig, ref in zip(counts_sig_ref[0], counts_sig_ref[1])] # calculate contrast (relative fluoresence intensity)
```

## Addition

### GUI software

A Graphical User Interface (GUI) version software based on this SDK is also implemented and has been used well in
our [QIM](https://quantum.lab.arizona.edu/) group. For more information please contact the author.

### Frequently Asked Questions

**Q: Where should I get to know more about Odmactor?**

**A:** If you have more demand or cooperation willingnes, please contact to
the [Quantum Information and Materials Group](https://quantum.lab.arizona.edu) of U Arizona.



### Copyright and License

Odmactor uses [The MIT license](LICENSE).

## Reference

1. https://www.degruyter.com/document/doi/10.1515/nanoph-2019-0209/html
2. https://zhuanlan.zhihu.com/p/361148655
