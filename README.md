# Ddmactor (ODMR Actor)
ODMR software SDK integrating functions of ODMR detection and spin manipulation

*This project is on continuous updating & using ...*

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
```python
# configure
channel_dict = {
    'laser': 1,
    'mw': 2,
    'apd': 3,
    'tagger': 5
}
tagger_input = {'apd': 1, 'asg': 2}

scheduler = PulseScheduler()  # Pulse-ODMR detection
scheduler.channel = channel_dict  # set ASG control channels
scheduler.tagger_input = tagger_input  # set tagger input channels
scheduler.configure_mw_paras(power=10)
scheduler.configure_odmr_seq(t_init, t_mw, t_read_sig=400, t_read_ref=400, inter_init_mw=inter_init_mw, N=N)
scheduler.set_mw_scan_freq_start_stop(freq_start, freq_end, freq_step)
scheduler.configure_tagger_counting(reader='cbm')  # 'Counter Between Markers' measurement

# run
scheduler.run()
```


