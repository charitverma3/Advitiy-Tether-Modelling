import time
import datetime
import os
import numpy as np

dir_now = datetime.datetime.now().strftime('%Y-%m-%d,%H:%M:%S')
dir_script = os.path.dirname(os.path.abspath(__file__))
dir_dest = os.path.join(dir_script, r'Logs', dir_now)
os.makedirs(dir_dest)

f_r = os.path.join(dir_dest, r'r.csv')

r = np.zeros((13,100))
np.savetxt(f_r,r,delimiter=',')