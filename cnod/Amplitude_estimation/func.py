import numpy as np
from qutip import*

vec = np.linspace(-2.5, 2.5, 51)

N = 40
a = tensor(qeye(2),destroy(N))
sz = tensor(sigmaz(),qeye(N))
sx = tensor(sigmax(),qeye(N))
sy = tensor(sigmay(),qeye(N))
sz = tensor(sigmaz(),qeye(N))

def D(beta):
    return ((beta*a.dag() - np.conjugate(beta)*a)).expm()


p0 = np.zeros((len(vec), len(vec)))
p1 = np.zeros((len(vec), len(vec)))

def probability(state, vec):
    for m, vx in enumerate(vec):
        for n, vy in enumerate(vec):
            state_d = D(vx + 1j*vy) * state
            p0[m, n] = abs(state_d.overlap(state))**2
            p1[m, n] = 1 - p0[m, n]        
    return p0, p1

def M_slope(data, vec, axis):
    "for the first and lase points:  p_(i+1)-p(i) ) / (alpha_(i+1)-alpha_(i))" 
    "the rest of points: p_(i+1)-p(i-1) ) / (alpha_(i+1)-alpha_(i-1))"
    step  = (vec.max() - vec.min()) / len(vec)
    if axis == 1:
        data_diff = data[:, 2:] - data[:, :-2]
        first_column = data[:,1] - data[:,0]
        last_column = data[:, -1] - data[:, -2]
        ## add first and last columns
        data_diff_add1 = np.insert(data_diff, 0, first_column, axis=1)
        data_diff_add2 = np.insert(data_diff_add1, data_diff_add1.shape[1], last_column, axis=1)
        slope = data_diff_add2 / (2 * step)
    else:
        data_diff = data[2:,:] - data[:-2, :]
        first_row = data[1, :] - data[0, :]
        last_row = data[-1, :] - data[-2, :]
        ## add first and last rows
        data_diff_add1 = np.insert(data_diff, 0, first_row, axis=0)
        data_diff_add2 = np.insert(data_diff_add1, data_diff_add1.shape[0], last_row, axis=0)
        slope = data_diff_add2 / (2 * step)
    return slope


def FI_calcluation(state, vec):
    xvec = vec
    yvec = vec
    p0 = probability(state, vec)[0]
    p1 = probability(state, vec)[1]
    f_11 = (M_slope(p0, xvec, axis=0) * M_slope(p0, xvec, axis=0)) / p0 + (M_slope(p1, xvec, axis=0) * M_slope(p1, xvec, axis=0)) / p1
    f_12 = (M_slope(p0, xvec, axis=0) * M_slope(p0, yvec, axis=1)) / p0 + (M_slope(p1, xvec, axis=0) * M_slope(p1, yvec, axis=1)) / p1
    f_22 = (M_slope(p0, yvec, axis=1) * M_slope(p0, yvec, axis=1)) / p0 + (M_slope(p1, yvec, axis=1) * M_slope(p1, yvec, axis=1)) / p1
    f_21 = f_12
    return f_11, f_12, f_21, f_22, f_11+f_22

