import numpy as np
import sys
from math import sqrt, floor
import matplotlib.pyplot as plt

class Simulation:
    def __init__(self, spin_flux, sim_params, medium_params, vacuum_params):
        self.spin_flux = spin_flux
        self.sim_p = sim_params
        self.jef_terms = self.sim_p["jef_terms"]
        self.dt = self.sim_p["dt"]
        self.med_p = medium_params
        self.vac_p = vacuum_params
        self.time = self.format_time()
        self.mediums = self.get_mediums()
        self.vacuum = self.get_vacuum()

    def format_time(self):
        time = np.arange(0, self.spin_flux.shape[0], self.dt)
        return time
    
    def get_mediums(self):
        N_time = len(self.time)
        mediums = []
        for medium in self.med_p:
            z_0, z_end = medium["z"]
            medium_spin_flux = self.spin_flux[:,z_0:z_end]
            mediums.append(Medium(medium, self.sim_p, medium_spin_flux))
        return mediums
    
    def get_vacuum(self):
        N_time = len(self.time)
        vac = Vacuum(self.vac_p, self.sim_p, N_time, self.mediums)
        return vac
                            
    def run(self, disp_progress = False):
        dt = self.sim_p["dt"]
        max_shift = floor(self.vacuum.min_max_delay[0][0])
        for t_i, t in enumerate(self.time):
            self.step(t_i, max_shift)
            if disp_progress:
                print(f"Step {t_i+1} of {len(self.time} completed")

    def step(self, t_i, max_shift):
        dt = self.sim_p["dt"]
        t = t_i * dt
        dV = self.sim_p["dx"] * self.sim_p["dy"] * self.sim_p["dz"]
        t_0 = t_i + max_shift - 2
        t_0 = 0 if t_0 < 0 else t_0
        t_len = t_i - t_0 + 1
        
        J_t_factor = 16022   # prefactor of Jefimenkos equations, converting from particle, nm, fs
                                    # to Coulomb, meter, second
        if self.jef_terms["rho"] or self.jef_terms["rho_t"]:
            use_dist = True

        for m_i, medium in enumerate(self.mediums):
            
            Nx, Ny, Nz = medium.N_points


            if medium.params["use_rho"]:
                pass
                # J_rho = ...


            if medium.params["use_Jx"]:
                J = np.zeros((t_len, Nx, Ny, Nz))
                
                # print(Nx, Ny, Nz)
                # print(t_len)
                exc_reg = np.tile(medium.exc_region, (t_len, 1, 1, 1))
                
                spin_flux = medium.spin_flux[t_0:t_i+1,:].reshape((t_len, 1, 1, Nz))
                # print(spin_flux)
                charge_flux = exc_reg * spin_flux * medium.params["theta"]
                # print(charge_flux.shape)
                for z_i, z in enumerate(self.vacuum.z):
                    #print(f"Distance {z_i + 1} of {len(self.vacuum.z) + 1} started")
                    Ex = 0
                    delta_t, delta_r = self.vacuum.time_delays[z_i][m_i], self.vacuum.source_vectors[z_i][m_i]
                    for index, val in np.ndenumerate(delta_t):
                        # print("----")
                        x_i, y_i, z_i2 = index
                        t_shift = delta_t[index]
                        r = delta_r[index]
                        t_i_ret = (t_len - 1 + t_shift)/dt
                        charge_flux_xyz = charge_flux[:,x_i, y_i, z_i2].copy().flatten()
                        # print(charge_flux_xyz)
                        # raise ValueError
                        J_t = get_J_t(t_i_ret, dt, charge_flux_xyz)
                        Ex += -J_t*J_t_factor*dV/r
                    self.vacuum.Ex_array[t_i, z_i] = Ex

    def save_Ex(self, path):
        arr = sim.vacuum.Ex_array.copy()
        np.savetxt(path, arr)

class Medium:
    def __init__(self, params, sim_params, spin_flux):
        self.params = params
        self.spin_flux = spin_flux
        self.format_domain(sim_params)
        self.N_points = self.get_N_points(sim_params)
        # self.arrays = self.get_arrays(N_time)
        self.exc_region = self.get_excitation_region(sim_params)
        self.use = self.get_usage()
    
    def format_domain(self, sim_params):
        for coord in ["x", "y", "z"]:
            old_min, old_max = self.params[coord]
            shift = sim_params[f"d{coord}"]/2
            self.params[coord] = (old_min + shift, old_max - shift)

    def get_N_points(self, sim_p):
        N_points = lambda coord: int((self.params[coord][1] - self.params[coord][0])/sim_p[f"d{coord}"] + 1 )
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
        for coord in ["x", "y", "z"]:
            if self.params[f"use_J{coord}"]:
                return True
        if self.params["use_rho"] or self.params["use_field"]:
            return True
        else:
            return False

class Vacuum:
    def __init__(self, params, sim_params, N_time, mediums):
        self.params = params
        self.z = self.get_z_points()
        self.Ex_array = self.get_array(N_time)
        self.time_delays, self.source_vectors = self.get_time_and_source(sim_params, mediums)
        self.min_max_delay = self.get_min_max_delay()
    
    def get_z_points(self):
        return self.params["z"]

    def get_array(self, N_time):
        arr = np.zeros((N_time, len(self.z)))
        return arr

    def get_time_and_source(self, sim_p, mediums):
        mult_factor = 10**6 / 299792458  # Division by C and conversion from ns to fs
        # print(mult_factor)
        all_time_delays = []
        all_source_vectors = []
        for z in self.z:
            delta_t_max = False
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
                    medium_delta_t_max = np.amax(delta_t)
                    if not delta_t_max:
                        delta_t_max = medium_delta_t_max
                    elif delta_t_max < medium_delta_t_max:
                        delta_t_max = medium_delta_t_max
            arr_type = type(np.zeros(1))
            shifted_time_delay_z = []
            for delta_t in time_delay_z:
                if isinstance(delta_t, arr_type):
                    # print("yo")
                    s = delta_t.shape
                    shift = np.ones(s) * delta_t_max
                    # print(shift)
                    # print(delta_t)
                    shifted_delta_t = delta_t - shift
                    # print(shifted_delta_t)
                    shifted_time_delay_z.append(shifted_delta_t)
                else:
                    shifted_time_delay_z.append(None)
            # print(shifted_time_delay_z)
            all_time_delays.append(shifted_time_delay_z)
            all_source_vectors.append(source_vector_z)
        return all_time_delays, all_source_vectors

    def get_min_max_delay(self):
        min_max_list = []
        arr_type = type(np.zeros(1))
        for time_delay_z in self.time_delays:
            no_value = True
            for medium_delay in time_delay_z:
                if isinstance(medium_delay, arr_type):
                    medium_min = np.amin(medium_delay)
                    medium_max = np.amax(medium_delay)
                    if no_value:
                        min, max = medium_min, medium_max
                        no_value = False
                    else:
                        if medium_min < min:
                            min = medium_min
                        if medium_max > max:
                            max = medium_max
            min_max_list.append((min, max))
        
        # print(min_max_list)
        return min_max_list

def get_J(t_i_ret, time_series):
    if t_i_ret <= 0 :
        return 0
    else:
        # x_i, y_i, z_i = pos_index
        # print(arr.shape)
        # time_series = arr[:,x_i,y_i,z_i].copy()
        # print(time_series.shape)
        # print(time_series)
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

if __name__ == "__main__":
    # np.set_printoptions(precision=0)
    
    from input_params import sim_params, medium_params, vacuum_params, spin_flux

    sim = Simulation(spin_flux, sim_params, medium_params, vacuum_params)
    sim.run(disp_progress=True)
    E0 = sim.vacuum.Ex_array[:,0]
    sim.save_Ex("testsave.out")
    time = sim.time
    plt.plot(time, E0)
    plt.savefig("images/onething.jpg")
    # print(sim.mediums[1].arrays["Jx"][100,:,:,:])
    # print(sim.J_s[100,:]*0.068)
    
    # # print(sim.mediums[1].arrays)


