import numpy as np
import os
import sys
# p = os.getcwd()
# p = os.path.dirname(__file__)
sys.path.insert(0, r'/Users/giorgiocanalella/Downloads/analyze/fit_funcs')
# sys.path.insert(0, r'/Users/giorgiocanalella/Desktop/2 - NUS/0 - Modules/4 - Y2S2//Users/giorgiocanalella/Desktop/2 - NUS/0 - Modules/4 - Y2S2/PC2288 - Basic UROPS in Physics/QCrew/Code/analyze/fit_funcs')
# sys.path.insert(0, r'/Users/giorgiocanalella/Downloads/analyze/fit_funcs')
# print(p)
import sine as sine


def func(xs, amp=-1, tau=1e5, f0=5e-5, phi=np.pi/4, ofs=0.8, tP=1e6):
    return (0.5 * amp * np.exp(-xs / tau) * (1 + np.exp(-xs / tP) * np.sin(2*np.pi*xs*f0 + phi) ) ) + ofs

def guess(xs, ys):
    d = sine.guess(xs, ys)
    d['tau'] = (np.average(xs), 0, 10*xs[-1])
    d['tP'] = (np.average(xs), 0, 10*xs[-1])
    d['f0'] = 4.4e-5
    return d
