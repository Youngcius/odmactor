# Optically Detected Magnetic Resonance

## Continuous-wave ODMR

### Necessary parameters

- ``

```python
channel_dict = {
    'laser': 1,
    'mw': 2,
    'apd': 3,
    'tagger': 5
}
tagger_input = {'apd': 1, 'asg': 2}
```

## Pulse ODMR

Wave form for single period:

```
asg laser channel:
-----               ---------
|   |               |       |
|   |---------------|       |----

asg microwave channel:
        -------------
        |           |
--------|           |------------

asg tagger acquisition channel:
                    ---   ---
                    | |   | |
--------------------| |---| |----
```

![](../../test/pulse-seq.png)


### Necessary parameters
Unit of the following time parameters is "ns".
- `t_init`: time span for laser initialization, e.g. 5000
- `t_mw`: time span for microwave actual operation in a ASG period, e.g. 800
- `t_read_sig`: time span for fluorescence signal readout, e.g. 400
- `t_read_ref`: time span for reference signal readout, e.g. 400 if the parameter is not assigned, means single-pulse readout
- `inter_init_mw`: time interval between laser initialization and MW operation pulses, e.g. 3000
- `inter_readout`: interval between single readout pulse and reference signal readout, e.g. 200 when t_read_ref is assigned, this parameter will play its role
- `inter_period`: interval between two neighbor periods, e.g. 200
- `N`: number of ASG operation periods





