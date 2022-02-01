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
    Plot counting contrast figure along with frequencies
    :param freqs: frequencies, unit: GHz
    :param contrast: contrast, range between 0 and 1
    :param fname: if assigned, a '.png' file will be saved
    """
    plt.style.use('seaborn')
    plt.figure()
    plt.plot(freqs, contrast, 'o-')
    plt.xlabel('frequency (GHz)')
    plt.ylabel('count contrast')
    plt.title('Contrast')
    if fname is not None:
        plt.savefig(fname, dpi=400)
    plt.show()
