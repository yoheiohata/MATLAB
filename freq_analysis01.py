# encoding: utf-8
import csv
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
import pandas as pd

# samplingrate = 4 kHz
def calc_fft(data, samplerate):
    spectrum = fftpack.fft(data)
    amp = np.sqrt((spectrum.real ** 2) + (spectrum.imag ** 2))
    amp = amp /(len(data) / 2)
    phase = np.arctan2(spectrum.imag, spectrum.real)
    phase =np.degrees(phase)
    freq = np.linspace(0, samplerate, len(data))
    return spectrum, amp, phase, freq


def main():
    # 時系列逐次処理：csvデータを一行読んで一行捨てるため，メモリを圧迫しない

    # csvファイルを開く
    with open('freq_analysis.csv', mode="r") as f:

        # ヘッダ行のスキップ
        skip_row_count = 3
        while skip_row_count != 0:
            row_data = f.readline()
            skip_row_count -= 1

        # データの読み込み開始
        # 最初の1行目はループ前に読み込む
        row_data = f.readline()
        split_data = row_data.split(",")

        # 2行目以降の読み出し処理を開始
        dt_cnt = 0
        while row_data:
            t = float(split_data[0])
            r = float(split_data[1])

            dt_cnt+=1
            row_data = f.readline()
            split_data = row_data.split(",")

    print(calc_fft(3, 4000))

    # 関数を実行してcsvファイルをフーリエ変換するだけの関数を実行
    # df, df_fft = csv_fft(in_file='freq_analysis.csv', out_file='fft.csv')
            

if __name__ == "__main__":
    main()