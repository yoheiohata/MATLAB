from scipy.fft import fft, fftfreq
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging
import datetime
import os


def FFT(y, samplingrate):
    N = len(y)
    y = np.array(y)
    yf = fft(y)
    amp = np.abs(yf[0:N//2])
    amp = amp/(N/2)
    phase = np.arctan2(yf.imag, yf.real)
    phase = np.degrees(phase)
    # freq = np.linspace(0, samplingrate, N)
    freq = fftfreq(N, samplingrate)[:N//2]
    result = [amp, phase, freq]
    return result

def main():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(
        filename=f'./log/{os.path.basename(__file__)}_{datetime.datetime.today()}.log'
        )
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)8s %(message)s"))
    logger.addHandler(handler)

    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    fig = plt.figure(figsize=(10, 5))
    fig.tight_layout()
    ax = fig.add_subplot(121)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    # ax.legend()


    """ parameter """
    sample_rate = 1/4000
    with open("freq_analysis.csv", "r") as f:

        skip_row_count = 3
        while skip_row_count != 0:
            raw_datas = f.readline()
            skip_row_count -= 1
        
        raw_datas = f.read().split("\n")

    t = []
    r = []
    for row_data in raw_datas:
        if len(row_data) == 0:
            continue
        row = row_data.split(",")
        t.append(row[0])
        r.append(row[1])
    
    r = np.array(r)
    amp_phase_freq = FFT(r, sample_rate)
    amp = amp_phase_freq[0]
    phase = amp_phase_freq[1]
    freq = amp_phase_freq[2]
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Amplitude [r/min]')
    plt.plot(freq, amp, label="FFT")
    plt.legend()
    # plt.savefig('./figure/fft.pdf')
    ax.set_ylim()
    plt.show()


if __name__ == "__main__":
    main()