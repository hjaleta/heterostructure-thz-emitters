import numpy as np
from math import sqrt

from numpy.lib.index_tricks import index_exp

class Simulation:
    def __init__(self, J_s, sim_params, medium_params, vacuum_params):
        self.c = {"c": 299792458, "eps0": 8.85418781*10**-12}
        self.J_s = J_s
        self.sim_p = sim_params
        self.jef_terms = self.sim_p["jef_terms"]
        self.dt = self.sim_p["dt"]
        self.med_p = medium_params
        self.vac_p = vacuum_params
        self.time = self.format_time()
        self.mediums = self.get_mediums()
        self.vacuum = self.get_vacuum()

    def format_time(self):
        time = np.arange(0, self.J_s.shape[0], self.dt)
        return time
    
    def get_mediums(self):
        N_time = len(self.time)
        mediums = []
        for medium in self.med_p:
            mediums.append(Medium(medium, self.sim_p, N_time))
        return mediums
    
    def get_vacuum(self):
        N_time = len(self.time)
        vac = Vacuum(self.vac_p, self.sim_p, N_time, self.mediums)
        return vac
                            
    def run(self):
        dt = self.sim_p["dt"]
        if self.jef_terms["J"] or self.jef_terms["J_t"]:
            z_0 = 0
            for medium in self.mediums:     # This loop calculates the charge current induced by the ISHE
                z_end = z_0 + medium.N_points[2]
                if medium.params["use_J"]["x"]:
                    for t_i in range(len(self.time)):
                        k = medium.exc_region.copy() * self.J_s[t_i,z_0:z_end] * medium.params["theta"]
                        medium.arrays["Jx"][t_i,:,:,:] =  k
                z_0 = z_end    

        for t_i, t in enumerate(self.time):
            self.step(t_i)

    def step(self, t_i):
        dt = self.sim_p["dt"]
        t = t_i * dt
        dV = self.sim_p["dx"] * self.sim_p["dy"] * self.sim_p["dz"]
        z_0 = 0
        J_t_factor = 16022   # prefactor of Jefimenkos equations, converting from particle, nm, fs
                                    # to Coulomb, meter, second 
        
        if self.jef_terms["rho"] or self.jef_terms["rho_t"]:
            use_dist = True

        for z_i, z_vacuum in self.vacuum.z:
            Ex = 0
            for m_i, medium in enumerate(self.mediums):
                if medium.use:
                    label = medium.params["label"]
                    delta_t = self.vacuum.time_delays[z_vacuum][m_i].copy()
                    delta_r = self.vacuum.source_vectors[z_vacuum][m_i].copy()
                    z_end = int((medium.params["z"][1] - z_0)/self.sim_p["dz"])
                    for arr in medium.arrays:
                        if isinstance(medium.arrays[arr], type(np.array([0]))):
                            if arr == "Jx":
                                print(medium.arrays[arr])
                                # raise ValueError()
                                for index, val in np.ndenumerate(medium.arrays[arr][t_i,:,:,:]):
                                    t_shift = delta_t[index]
                                    t_i_ret = (t + t_shift)/dt
                                    print(t_i_ret)
                                    J_t = get_J_t(t_i_ret, dt, index, medium.arrays[arr])
                                    r = delta_r[index]
                                    Ex += J_t * J_t_factor / r
                            elif arr == "Jy":
                                pass
                            elif arr == "Jz":
                                pass
                            elif arr == "rho":
                                pass
                        

                                        

            
            if self.jef_terms["J"]:
                pass
            if self.jef_terms["J_t"]:
                pass

class Medium:
    def __init__(self, params, sim_params, N_time):
        self.params = params
        self.format_domain(sim_params)
        self.N_points = self.get_N_points(sim_params)
        self.arrays = self.get_arrays(N_time)
        self.exc_region = self.get_excitation_region(sim_params)
        self.use = self.get_usage()
    
    def format_domain(self, sim_params):
        for coord in ["x", "y", "z"]:
            old_min, old_max = self.params[coord]
            shift = sim_params[f"d{coord}"]/2
            self.params[coord] = (old_min + shift, old_max - shift)

    def get_N_points(self, sim_p):
        N_points = lambda coord: int((self.params[coord][1] - self.params[coord][0])/sim_p[f"d{coord}"] + 1 )# - int(self.use_FEM))
        Nx = N_points("x")
        Ny = N_points("y")
        Nz = N_points("z")
        return (Nx, Ny, Nz)

    def get_arrays(self, N_time):
        Nx, Ny, Nz = self.N_points
        arrays = {}
        
        for coord in self.params["use_J"]:
            if self.params["use_J"][coord]:
                arrays[f"J{coord}"] = np.zeros((N_time, Nx,Ny,Nz))
            else:
                arrays[f"J{coord}"] = None
        if self.params["use_rho"]:
            arrays["rho"] = np.zeros((N_time, Nx,Ny,Nz))
        else:
            arrays["rho"] = None
        for coord in ["x", "y", "z"]:
            if self.params["use_field"]:
                if self.params["use_field"][f"E{coord}"]:
                    arrays[f"E{coord}"] = np.zeros((N_time, Nx,Ny,Nz))
                else:
                    arrays[f"E{coord}"] = None
                if self.params["use_field"][f"B{coord}"]:
                    arrays[f"B{coord}"] = np.zeros((N_time, Nx,Ny,Nz))
                else:
                    arrays[f"B{coord}"] = None
        
        return arrays

    def get_excitation_region(self, sim_p):
        Nx, Ny, Nz = self.N_points
        exc_region = sim_p["exc_region"]
        
        x_min, x_max = self.params["x"]
        y_min, y_max = self.params["y"]
        dx = sim_p["dx"]
        dy = sim_p["dy"]
        X = np.arange(x_min, x_max + dx, dx)
        Y = np.arange(y_min, y_max + dy, dy)
        if (len(X) != Nx) or (len(Y) != Ny): 
            raise ValueError("Dimensions of excitation region is inconsistent")
        if exc_region["shape"] == "rectangle":
            x_side, y_side = exc_region["x_side"]/2, exc_region["y_side"]/2
            l1 = []
            for x_i in X:
                l2 = []
                for y_i in Y:
                    if (abs(x_i) <= x_side) and (abs(y_i) <= y_side):
                        l2.append(1)
                    else:
                        l2.append(0)
                l1.append(l2)

        elif exc_region["shape"] == "circle":
            r_squared = exc_region["radius"]**2
            l1 = []
            for x_i in X:
                l2 = []
                for y_i in Y:
                    if x_i**2 + y_i**2 < r_squared:
                        l2.append(1)
                    else:
                        l2.append(0)
                l1.append(l2)
        else:
            raise ValueError("Invalid shape of excitation region, use 'circle' or 'square'")
        region = np.array([l1]*Nz)
        region = np.moveaxis(region, 0, 2)
        # print(region.shape)
        return region

    def get_usage(self):
        # self.use = boo
        for coord in self.params["use_J"]:
            if self.params["use_J"][coord]:
                return True
        if self.params["use_rho"] or self.params["use_field"]:
            return True
        else:
            return False

class Vacuum:
    def __init__(self, params, sim_params, N_time, mediums):
        self.params = params
        self.z = self.get_z_points()
        self.array = self.get_arrays(N_time)
        self.time_delays, self.source_vectors = self.get_time_and_source(sim_params, mediums)
    
    def get_z_points(self):
        return self.params["z"]

    def get_arrays(self, N_time):
        arr = np.zeros((N_time, len(self.z)))
        return arr

    def get_time_and_source(self, sim_p, mediums):
        mult_factor = 10**6 / 299792458  # Division by C and conversion from ns to fs
        all_time_delays = []
        all_source_vectors = []
        for z in self.z:
            time_delay_z =  [] 
            source_vector_z = []
            for medium in mediums:
                if not medium.use:
                    time_delay_z.append(None)
                    source_vector_z.append(None)
                else:
                    x_min, x_max = medium.params["x"]
                    y_min, y_max = medium.params["y"]
                    z_min, z_max = medium.params["z"]
                    delta_t = np.zeros(medium.N_points)
                    delta_r = np.zeros(medium.N_points)
                    dx, dy, dz =  sim_p["dx"], sim_p["dy"], sim_p["dz"]
                    for index, val in np.ndenumerate(delta_t):
                        delta_x = x_min + index[0] * dx
                        delta_y = y_min + index[1] * dy
                        delta_z = z - (z_min + index[2] * dz) 
                        delta_r[index] = sqrt(delta_x**2 + delta_y**2 + delta_z**2)
                        # source_vector[index] = delta_r
                    
                    delta_t = - delta_r * mult_factor
                        
                    time_delay_z.append(delta_t)
                    source_vector_z.append(delta_r)
            all_time_delays.append(time_delay_z)
            all_source_vectors.append(source_vector_z)
        return time_delays, source_vectors

def get_J(t_i_ret, pos_index, arr):
    if t_i_ret <= 0:
        return 0
    else:
        x_i, y_i, z_i = pos_index
        time_series = arr[:,x_i,y_i,z_i].copy()
        num = t_i_ret % 1
        t_0, t_1 = int(round(t_i_ret - num)), int(round(t_i_ret-num + 1)), 
        J = interpolation(time_series[t_0], time_series[t_1], num)
        return J

def get_J_t(t_i_ret, dt, pos_index, arr):
    J_2 = get_J(t_i_ret, pos_index, arr)
    J_1 = get_J(t_i_ret - 1, pos_index, arr)
    J_0 = get_J(t_i_ret - 2, pos_index, arr)
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
    


if __name__ == "__main__":
    np.set_printoptions(precision=0)
    
    from input_params import  sim_params, medium_params, vacuum_params, spin_current

    sim = Simulation(spin_current, sim_params, medium_params, vacuum_params)
    sim.run()
    print(sim.mediums[1].arrays["Jx"][100,:,:,:])
    print(sim.J_s[100,:]*0.068)
    
    # print(sim.mediums[1].arrays)


