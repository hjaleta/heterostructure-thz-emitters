import numpy as np
from math import sqrt

def get_J(t_i_ret, time_series):
    if t_i_ret <= 0 :
        return 0
    else:
        if np.abs(t_i_ret - round(t_i_ret)) < 0.001:
            J = time_series[round(t_i_ret)]
        else:
            num = t_i_ret % 1
            t_0, t_1 = round(t_i_ret - num), round(t_i_ret-num + 1)
            if len(time_series) > 1:
                J = interpolation(time_series[t_0], time_series[t_1], num)
            else:
                J = time_series[0]
        return J

def get_r_tuple(x_s,y_s,z_s, z_target, interfaces):
    r = []
    for interface in interfaces:
        if z_s < interface:
            z_i = interface
            factor = 1 - (z_i-z_s)/(z_target-z_s)
            # print(factor)
            y_i = y_s*factor
            x_i = x_s*factor
            delta_r = sqrt((z_i-z_s)**2 + (x_i-x_s)**2 + (y_i-y_s)**2)
            r.append(delta_r)
            z_s = z_i
            y_s = y_i
            x_s = x_i
            # print(x_s, y_s, z_s)
        else:
            r.append(0)
    delta_r = sqrt((z_s-z_target)**2 + x_s**2 + y_s**2)
    r.append(delta_r)
    return np.array(r)

def calc_delta_t(r_tuple, n_real):
    r_effective = 0
    for r, n in zip(r_tuple, n_real):
        r_effective += r*n
    t = r_effective/299 # Speed of light in [nm / fs]

    return t

def get_J_t(t_i_ret, dt, time_series):
    J_2 = get_J(t_i_ret, time_series)
    J_1 = get_J(t_i_ret - 1, time_series)                                                                                                                                                 
    J_0 = get_J(t_i_ret - 2, time_series)
    J_t = back_diff2(J_2, J_1, J_0, dt)
    return J_t

def interpolation(val1, val2, number):
    if (number < 0) or (number > 1):
        raise ValueError("number must be in interval [0,1]")
    num = val1 + number*(val1 - val2)
    return num

def back_diff1(y1, y0, dt):
    yprim = (y1 -y0)/dt
    return yprim

def back_diff2(y2, y1, y0, dt):
    yprim = (1.5*y2 - 2*y1 + 0.5*y0)/dt
    return yprim
