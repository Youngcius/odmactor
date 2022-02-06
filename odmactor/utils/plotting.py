"""
Utils plotting functions
"""

from time import time
from typing import List
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def plot_pulse_odmr():
    pass


# TODO: 合理的公式拟合？？有没有必要呢

def plot_ramsey(times: List[float], contrast: List[float], fname: str = None) -> Figure:
    """
    Plot Ramsey detection figure along with free precession time intervals
    :param times: frequencies, unit: ns
    :param contrast: contrast, range between 0 and 1
    :param fname: if assigned, a '.png' file will be saved
    :return: a matplotlib Figure instance
    """
    fig = plt.figure()
    plt.plot(times, contrast, 'o-')
    plt.title('Ramsey detection for T2*')
    plt.xlabel('Free precession time (ns)')
    plt.ylabel('Contrast $N_{sig}/N_{ref}$')
    if fname is not None:
        plt.savefig(fname, dpi=350)
    return fig


def plot_rabi(times: List[float], contrast: List[float], fname: str = None) -> Figure:
    """
    Plot Rabi assilication figure along with MW pulse duration
    :param times: frequencies, unit: ns
    :param contrast: contrast, range between 0 and 1
    :param fname: if assigned, a '.png' file will be saved
    :return: a matplotlib Figure instance
    """
    fig = plt.figure()
    plt.plot(times, contrast, 'o-')
    plt.title('Rabi oscillation for calibrated $\pi$ pulse')
    plt.xlabel('MW pulse duration (ns)')
    plt.ylabel('Contrast $N_{sig}/N_{ref}$')
    if fname is not None:
        plt.savefig(fname, dpi=350)
    return fig


def plot_t1(times: List[float], contrast: List[float], fname: str = None) -> Figure:
    """
    Plot T1 relaxation figure along with laser delay time intervals
    :param times: frequencies, unit: ns
    :param contrast: contrast, range between 0 and 1
    :param fname: if assigned, a '.png' file will be saved
    :return: a matplotlib Figure instance
    """
    fig = plt.figure()
    plt.plot(times, contrast, 'o-')
    plt.title('T1 Relaxometry')
    plt.xlabel('Relaxation time (ns)')
    plt.ylabel('Contrast $N_{sig}/N_{ref}$')
    if fname is not None:
        plt.savefig(fname, dpi=350)
    return fig


def plot_freq_contrast(freqs, contrast, fname: str = None) -> Figure:
    """
    Plot counting contrast figure along with frequencies
    :param freqs: frequencies, unit: GHz
    :param contrast: contrast, range between 0 and 1
    :param fname: if assigned, a '.png' file will be saved
    :return: a matplotlib Figure instance
    """
    plt.style.use('seaborn')
    fig = plt.figure()
    plt.plot(freqs, contrast, 'o-')
    plt.xlabel('frequency (GHz)')
    plt.ylabel('count contrast')
    plt.title('Contrast')
    if fname is not None:
        plt.savefig(fname, dpi=350)
    return fig
