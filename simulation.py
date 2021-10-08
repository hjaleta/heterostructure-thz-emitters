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
        self.interfaces = self.get_interfaces()
        self.mediums = self.get_mediums()
        self.vacuum = self.get_vacuum()

    def format_time(self):
        time = np.arange(0, self.spin_flux.shape[0], self.dt)
        return time
    
    def get_interfaces(self):
        interfaces = [p["z"][1] for p in self.med_p]
        return interfaces

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
        vac = Vacuum(self.vac_p, self.sim_p, N_time, self.mediums, self.interfaces)
        return vac
                            
    def run(self, disp_progress = False):
        dt = self.sim_p["dt"]
        max_shift = floor(self.vacuum.min_max_delay[0][0])
        for t_i, t in enumerate(self.time):
            self.step(t_i, max_shift)
            if disp_progress:
                print(f"Step {t_i+1} of {len(self.time)} completed")

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
        arr = self.vacuum.Ex_array.copy()
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

class E_field:
    def __init__(self, z, sim_params, mediums, N_time, interfaces):
        self.z = z
        self.Ex = np.zeros(N_time)
        self.delta_r = self.get_delta_r(sim_params, mediums, interfaces)
        self.delta_t = self.get_delta_t(mediums)

    def get_delta_r(self, sim_params, mediums, interfaces):
        all_delta_r = []
        for m_i, medium in enumerate(mediums):
            if not medium.use:
                all_delta_r.append(None)
            else:
                x_min, x_max = medium.params["x"]
                y_min, y_max = medium.params["y"]
                z_min, z_max = medium.params["z"]
                #delta_t = np.zeros(medium.N_points)
                Nx, Ny, Nz = medium.N_points
                delta_r = np.zeros((Nx,Ny,Nz, len(mediums)+1))
                dx, dy, dz =  sim_params["dx"], sim_params["dy"], sim_params["dz"]
                for (x_i,y_i,z_i), val in np.ndenumerate(delta_r[:,:,:,0]):
                    x_s = x_min + x_i * dx
                    y_s = y_min + y_i * dy
                    z_s = z_min + z_i * dz
                    r = get_r_tuple(x_s, y_s, z_s, self.z, interfaces)
                    delta_r[x_i,y_i,z_i,:] = r
                all_delta_r.append(delta_r)
        return all_delta_r

    def get_delta_t(self, mediums):
        n = []
        for medium in mediums:
            n.append(medium.params["n"])
        n.append(1)
        all_delta_t = []
        for arr in self.delta_r:
            if arr is None:
                all_delta_t.append(None)
            else:
                Nx,Ny,Nz,Nm = arr.shape
                delta_t = np.zeros((Nx,Ny,Nz))
                for (x_i, y_i, z_i), _ in np.ndenumerate(delta_t):
                    tup = arr[x_i,y_i,z_i,:].copy()
                    t = calc_delta_t(tup, n)
                    delta_t[x_i,y_i,z_i] = t
                t_min = np.amin(delta_t)
                delta_t -= np.ones(delta_t.shape)*t_min
                all_delta_t.append(delta_t)
        return all_delta_t
         
class Vacuum:
    def __init__(self, params, sim_params, N_time, mediums, interfaces):
        self.params = params
        self.z = self.get_z_points()
        self.E_fields = self.get_E_fields(sim_params, mediums, N_time, interfaces)

    def get_z_points(self):
        return self.params["z"]

    def get_E_fields(self, sim_params, mediums, N_time, interfaces):
        E_fields = []
        for z in self.z:
            E = E_field(z, sim_params, mediums, N_time, interfaces)
            E_fields.append(E)
        return E_fields

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
    # factor = 1/(3 * 10**2)
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

if __name__ == "__main__":
    # np.set_printoptions(precision=0)
    
    from input_params import sim_params, medium_params, vacuum_params, spin_flux

    sim = Simulation(spin_flux, sim_params, medium_params, vacuum_params)
    print(sim.vacuum.E_fields[0].delta_r)
    print(sim.vacuum.E_fields[0].delta_t[1])
    print(sim.vacuum.E_fields[0].delta_t[1].shape)
    # sim.run(disp_progress=True)
    # E0 = sim.vacuum.Ex_array[:,0]
    # sim.save_Ex("testsave.out")
    # time = sim.time
    # plt.plot(time, E0)
    # plt.savefig("images/onething.jpg")
    # print(sim.mediums[1].arrays["Jx"][100,:,:,:])
    # print(sim.J_s[100,:]*0.068)
    
    # # print(sim.mediums[1].arrays)

