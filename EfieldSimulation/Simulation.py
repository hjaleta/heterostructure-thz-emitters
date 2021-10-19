import numpy as np
import sys
from time import perf_counter
from math import sqrt, floor, ceil
import matplotlib.pyplot as plt
from E_fieldSimulation.help_functions import (get_J, get_J_t, get_r_tuple, 
            back_diff1, back_diff2, interpolation, calc_delta_t)

class Simulation:
    def __init__(self, spin_flux, sim_params, medium_params, vacuum_params, name = ""):
        self.name = name
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
        mediums = []
        for medium in self.med_p:
            # print(medium)
            z_0, z_end = medium["z"]
            # print(z_0,z_end)
            medium_spin_flux = self.spin_flux[:,z_0:z_end]
            mediums.append(Medium(medium, self.sim_p, medium_spin_flux))
        return mediums

    def get_vacuum(self):
        N_time = len(self.time)
        vac = Vacuum(self.vac_p, self.sim_p, N_time, self.mediums, self.interfaces)
        return vac
                            
    def run(self, print_percent = False):
        first_time = perf_counter()
        
        dt = self.sim_p["dt"]
        if print_percent:
            percent_step = ceil((len(self.time)*print_percent)/100)
        percent = -print_percent
        max_shift = ceil(self.vacuum.max_delta_t)
        for t_i, t in enumerate(self.time):
            
            self.step(t_i, max_shift)
            if print_percent and t_i % percent_step == 0:
                percent += print_percent
                print(f"Simulation {self.name}\n{percent} % completed")
                current_time = perf_counter()
                print(f"Total simulation time: {current_time - first_time:.0f} s\n----------")
        current_time = perf_counter()
        print(f"Simulation {self.name} completed in {current_time-first_time:.0f} s\n----------")

    def step(self, t_i, max_shift):
        dt = self.sim_p["dt"]
        t = t_i * dt
        dV = self.sim_p["dx"] * self.sim_p["dy"] * self.sim_p["dz"]
        t_0 = t_i - max_shift - 2
        t_0 = 0 if t_0 < 0 else t_0
        t_len = t_i - t_0 + 1
        
        J_t_factor = 16022   # prefactor of Jefimenkos equations, converting from particle, nm, fs
                                    # to Coulomb, meter, second

        for m_i, medium in enumerate(self.mediums):
            
            Nx, Ny, Nz = medium.N_points


            if medium.params["use_rho"]:
                pass

            if medium.params["use_Jx"]:
                # print("WOO")
                # J = np.zeros((t_len, Nx, Ny, Nz))
                exc_reg = np.tile(medium.exc_region, (t_len, 1, 1, 1))
                spin_flux = medium.spin_flux[t_0:t_i+1,:].reshape((t_len, 1, 1, Nz))
                charge_flux = exc_reg * spin_flux * medium.params["theta"]
                # print(np.sum(charge_flux, axis = 1))
                for E in self.vacuum.E_fields:
                    if E.delta_r[m_i] is None:
                        continue
                    Ex = 0
                    delta_t, delta_r_tot = E.delta_t[m_i], E.delta_r_tot[m_i]
                    # z_target = E.z
                    for index, t_shift in np.ndenumerate(delta_t):
                        x_i, y_i, z_i = index
                        r = delta_r_tot[index]
                        t_i_ret = (t_len - 1 - t_shift)/dt
                        charge_flux_xyz = charge_flux[:,x_i, y_i, z_i].copy().flatten()
                        J_t = get_J_t(t_i_ret, dt, charge_flux_xyz)
                        # print(J_t)
                        Ex += -J_t*J_t_factor*dV/r
                    E.Ex[t_i] = Ex
                    # print(Ex)

    def save_Ex(self, path):
        all_arrs = []
        for E in self.vacuum.E_fields:
            E_data = E.Ex.copy().reshape((1,-1))
            all_arrs.append(E_data)
        arr = np.array(all_arrs)
        np.savetxt(path, arr)

class Medium:
    def __init__(self, params, sim_params, spin_flux):
        self.params = params
        self.spin_flux = spin_flux
        self.format_domain(sim_params)
        self.N_points = self.get_N_points(sim_params)
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
    def __init__(self, params, sim_params, N_time, mediums, interfaces):
        self.params = params
        self.z = self.get_z_points()
        self.E_fields = self.get_E_fields(sim_params, mediums, N_time, interfaces)
        self.max_delta_t = self.get_max_delta_t()
        
    def get_z_points(self):
        return self.params["z"]

    def get_E_fields(self, sim_params, mediums, N_time, interfaces):
        E_fields = []
        
        for z in self.z:
            E = E_field(z, sim_params, mediums, N_time, interfaces)
            E_fields.append(E)
        return E_fields

    def get_max_delta_t(self):
        maxes =[]
        for E in self.E_fields:
            maxes.append(E.delta_t_max)
        return max(maxes)

class E_field:
    def __init__(self, z, sim_params, mediums, N_time, interfaces):
        self.z = z
        self.Ex = np.zeros(N_time)
        self.delta_r = self.get_delta_r(sim_params, mediums, interfaces)
        self.delta_r_tot = self.get_delta_r_tot()
        self.delta_t, self.delta_t_max = self.get_delta_t(mediums)

    def get_delta_r(self, sim_params, mediums, interfaces):
        all_delta_r = []
        for m_i, medium in enumerate(mediums):
            if medium.use:
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
            else:
                all_delta_r.append(None)
        return all_delta_r

    def get_delta_r_tot(self):
        all_delta_r_tot = []
        for delta_r in self.delta_r:
            if delta_r is None:
                all_delta_r_tot.append(None)
            else:
                r_array = np.zeros(delta_r[:,:,:,0].shape)
                for (xi,yi,zi), _ in np.ndenumerate(delta_r[:,:,:,0]):
                    r_tot = np.sum(delta_r[xi,yi,zi,:])
                    r_array[xi,yi,zi] = r_tot
                all_delta_r_tot.append(r_array)
        return all_delta_r_tot

    def get_delta_t(self, mediums):
        n = []
        for medium in mediums:
            n.append(medium.params["n"])
        n.append(1)
        all_delta_t = []
        delta_t_maxes = []
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
                delta_t_maxes.append(np.amax(delta_t))
        delta_t_max = max(delta_t_maxes)

        return all_delta_t, delta_t_max

if __name__ == "__main__":
    pass
    
    # from freq_analysis.broadband import Signal
    # signal_params = {"padfactor":9, "n_interpol":1, "dt":1, "bb_fraction": 10}

    # sim = Simulation(spin_flux, sim_params, medium_params, vacuum_params)
    # sim_signals = {}

    # sim.run(disp_progress=True)
    # for E_i, E in enumerate(sim.vacuum.E_fields):
    #     key = str(E.z)
    #     sim_signals[key] = Signal(E.Ex, signal_params)
    # for z in sim_signals:
    #     print(sim_signals[z].bandwidth)
    #     sim_signals[z].plot_bb("nice_bb")
    
