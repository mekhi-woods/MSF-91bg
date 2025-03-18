import os
import glob
import utils  # Import of utils.py
import shutil
import datetime
import numpy as np
import time as systime
from astropy.table import Table, vstack


if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')