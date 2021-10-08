from math import sqrt
import numpy as np

def get_r_tuple(x_s,y_s,z_s, z_target, absorption_regions):
    r = []
    for region in absorption_regions:
        if z_s < region[1]:
            z_i = region[1]
            factor = -(z_i-z_s)/(z_target-z_s)
            y_i = y_s*factor
            x_i = x_s*factor
            delta_r = sqrt((z_i-z_s)**2 + (x_i-x_s)**2 + (y_i-y_s)**2)
            r.append(delta_r)
            z_s = z_i
            y_s = y_i
            x_s = x_i
        else:
            r.append(None)
    delta_r = sqrt((z_s-z_target)**2 + x_s**2 + y_s**2)
    r.append(delta_r)
    return np.array(delta_r)

x_s, y_s, z_s = 10,10,10
z_target = 30
absorption_regions = [(0,5,4), (5,12,3), (12,17, 4)]
r = get_r_tuple(x_s, y_s, z_s, z_target, absorption_regions)
print(r)