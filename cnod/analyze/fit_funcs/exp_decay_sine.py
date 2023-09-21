import numpy as np
import os
import sys
# sys.path.insert(0, r'/Users/giorgiocanalella/Desktop/2 - NUS/0 - Modules/4 - Y2S2//Users/giorgiocanalella/Desktop/2 - NUS/0 - Modules/4 - Y2S2/PC2288 - Basic UROPS in Physics/QCrew/Code/analyze/fit_funcs')
sys.path.insert(0, r'/Users/giorgiocanalella/Downloads/analyze/fit_funcs')
import sine as sine


def func(xs, amp=1, f0=0.05, phi=np.pi/4, ofs=0, tau=0.5):
    return amp * np.sin(2*np.pi*xs*f0 + phi) * np.exp(-xs / tau) + ofs

def guess(xs, ys):
    d = sine.guess(xs, ys)
    d['tau'] = (np.average(xs), 0, 10*xs[-1])
    return d
