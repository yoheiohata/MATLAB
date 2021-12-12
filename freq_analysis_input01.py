import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import datetime
import os

def log():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(
        filename = f'./log/{os.path.basename(__file__)}{datetime.datetime.today()}.log',
        mode='w'
        )
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)8s %(message)s"))
    logger.addHandler(handler)

def pltParams():
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

    pltParams()

    input_csvname = './freq_analysis.csv'
    input_csv = pd.read_csv(input_csvname, header=2)

    first_column_data = input_csv[input_csv.keys()[0]]
    second_column_data = input_csv[input_csv.keys()[1]]

    plt.xlabel(input_csv.keys()[0])
    plt.ylabel(input_csv.keys()[1])

    plt.plot(first_column_data, second_column_data)
    plt.show()

    log()



if __name__ == "__main__":
    main()