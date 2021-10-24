import matplotlib.pyplot as plt

"""
Utils plotting functions
"""


def plot_pulse_odmr():
    pass


def plot_ramsey():
    pass


def plot_freq_contrast(freqs, contrast, fname: str = None):
    """

    :param freqs: frequencies, unit: GHz
    :param contrast: contrast, range between 0 and 1
    """
    plt.style.use('seaborn')
    plt.plot(freqs, contrast)
    plt.xlabel('frequency (GHz)')
    plt.ylabel('count contrast')
    if fname is not None:
        plt.savefig(fname, dpi=400)
    plt.show()
