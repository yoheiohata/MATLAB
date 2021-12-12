from scipy.fft import fft, fftfreq
import numpy as np
import pandas as pd
import matplotlib
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
    result = [freq, amp, phase]
    return result

def log():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(
        filename=f'./log/{os.path.basename(__file__)}_{datetime.datetime.today()}.log'
        )
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)8s %(message)s"))
    logger.addHandler(handler)

def plt_params():
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    fig = plt.figure(figsize=(10, 5))
    fig.tight_layout()
    ax = fig.add_subplot(121)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

def main():
    plt_params()
    sample_rate = 1/4000

    input_csvname = './freq_analysis.csv'
    input_csv = pd.read_csv(input_csvname, header=2)

    first_column_data = input_csv[input_csv.keys()[0]]
    second_column_data = input_csv[input_csv.keys()[1]]

    first_column_data = np.array(first_column_data)
    second_column_data = np.array(second_column_data)

    freq_amp_phase = FFT(second_column_data, sample_rate)

    pd.DataFrame(freq_amp_phase).T.to_csv('./fft.csv', columns=[0,1] ,index=False)

    log()

if __name__ == "__main__":
    main()