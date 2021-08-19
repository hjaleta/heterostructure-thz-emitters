import numpy as np
from math import sqrt

class Simulation:
    def __init__(self, J_s, sim_params, medium_params, vacuum_params):
        self.c = {"c": 299792458}
        self.J_s = J_s
        self.sim_p = sim_params
        self.jef_terms = self.sim_p["jef_terms"]
        self.dt = self.sim_p["dt"]
        self.med_p = medium_params
        self.vac_p = vacuum_params
        self.time = self.format_time()
        self.mediums = self.get_mediums()
        self.vacuum = self.get_vacuum()
        self.exc_region = self.get_exc_region()
        self.get_time_delay()


    def format_time(self):
        time = np.arange(0, self.J_s.shape[0], self.dt)
        return time
    
    def get_mediums(self):
        N_time = len(self.time)
        mediums = []
        for medium in self.med_p:
            mediums.append(Medium(medium, self.sim_p, N_time))
        return mediums
    
    def get_exc_region(self):
        N_points = self.mediums[0].N_points
        N_z = self.J_s.shape[1]
        x_min, x_max = self.mediums[0].params["x"]
        y_min, y_max = self.mediums[0].params["y"]
        dx = self.sim_p["dx"]
        dy = self.sim_p["dy"]
        # print(x_min, x_max, dx)
        X = np.arange(x_min, x_max + dx, dx)
        Y = np.arange(y_min, y_max + dy, dy)
        if len(X) != N_points[0] or len(Y) != N_points[1]: 
            print(len(X), N_points[0])
            print(len(Y), N_points[1])
            raise ValueError("Bajs")

        if self.sim_p["exc_region"]["shape"] == "circle":
            r_squared = self.sim_p["exc_region"]["radius"]**2
            l1 = []
            for x_i in X:
                l2 = []
                for y_i in Y:
                    if x_i**2 + y_i**2 < r_squared:
                        l2.append(1)
                    else:
                        l2.append(0)
                l1.append(l2)
            
            region = np.array([l1]*N_z)
            region = np.moveaxis(region, 0, 2)
            # print(region.shape)
            return region

    def get_time_delay(self):
        mult_factor = 10**6 / self.c["c"]   # Division by C and conversion from ns to fs
        for medium1 in self.mediums:
            if medium1.params["use_field"]:
                time_delays = {}
                for medium2 in self.mediums:
                    if medium2.params["label"][0:3] in ["non", "Mag"]:
                        x_min, x_max = medium2.params["x"]
                        y_min, y_max = medium2.params["y"]
                        z_min, z_max = medium2.params["z"]
                        time_delay = np.zeros(medium2.N_points)
                        dx, dy, dz =  self.sim_p["dx"], self.sim_p["dy"], self.sim_p["dz"]
                        for z in np.arange(medium1.params["z"][0], medium1.params["z"][0], self.sim_p["dz2"]):
                            time_delay = medium2.array[0,:,:,:].copy()
                            for index, val in np.ndenumerate(time_delay):
                                delta_x = x_min + index[0] * dx
                                delta_y = y_min + index[1] * dy
                                delta_z = z - (z_min + index[2] * dz) 
                                delta_r = sqrt(delta_x**2 + delta_y**2 + delta_z**2)
                                delta_t = delta_r * mult_factor
                                time_delay[index] = delta_t
                            

                

    def run(self):
        dt = self.sim_p["dt"]
        if self.jef_terms["J"] or self.jef_terms["J_t"]:
            z_0 = 0
            for medium in self.mediums:     # This loop calculates the charge current induced by the ISHE
                z_end = z_0 + medium.N_points[2]
                if medium.params["use_J"]["x"]:
                    print(medium.arrays["Jx"][1,:,:,:].shape)
                    for t_i in range(len(self.time)):
                        k = self.exc_region.copy()[:,:,z_0:z_end] * self.J_s[t_i,z_0:z_end] * medium.params["gamma"]
                        medium.arrays["Jx"][t_i,:,:,:] =  k
                z_0 = z_end    


    def step(self, t_i):
        z_0 = 0
        for medium in self.mediums:
            z_end = int((medium.params["z"][1] - z_0)/self.sim_p["dz"])
            for arr in medium.arrays:
                if medium.arrays[arr] == None:
                    continue
                elif arr == "J":
                    medium.arrays[arr][t_i,:,:,:] = self.exc_region.copy() * self.J_s.copy()[t_i,z_0:z_end]
            if self.jef_terms["rho"]:
                pass
            if self.jef_terms["rho_t"]:
                pass
            if self.jef_terms["J"]:
                pass
            if self.jef_terms["J_t"]:
                pass

class Medium:
    def __init__(self, params, sim_params, N_time):
        self.params = params
        self.use_FEM = self.params['FEM']
        self.format_domain(sim_params)
        self.N_points = self.get_N_points(sim_params)
        self.arrays = self.get_arrays(N_time)
    
    def format_domain(self, sim_params):
        if self.use_FEM:
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
        if self.params["label"][0:3] in ["non", "Mag"]:
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
        elif self.params["label"][0:6] == "vacuum":
            for coord in self.params["use_J"]:
                arrays[f"J{coord}"] = None
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

class Vacuum:
    def __init__(self, params, sim_params, N_time):
        self.params = params
        self.N_points = self.get_N_points(sim_params)
        self.arrays = self.get_arrays(N_time)

if __name__ == "__main__":
    jef_terms = {"rho": False, "rho_t": False, "J": False, "J_t": True} # Which terms from Jefimenko's equations to include
    exc_region = {"shape": "circle", "radius": 1000} # The region where the laser pulse hits the plates
    sim_params = {"dx":500, "dy":500, "dz":1, "dz2": 1000, "dt": 1, # time is fs, length is nm
                "jef_terms": jef_terms, "exc_region": exc_region} 
    med_params = [
        {"label": "Mag",
        "x": (-2000, 2000), "y":(-2000, 2000), "z": (0,5), "FEM":True,
        "use_J":{"x": False, "y": False, "z":False}, 
        "use_rho":False, 
        "use_field": False,
        "material":"Fe", "gamma":0},

        {"label": "nonMag",
        "x": (-2000, 2000), "y":(-2000, 2000), "z": (5,8), "FEM":True,
        "use_J":{"x": True, "y": False, "z":False}, 
        "use_rho":False, 
        "use_field": False,
        "material":"Pt", "gamma":0.068},

        {"label":"vacuum",
        "x": (0, 0), "y":(0, 0), "z": (10, 10010), "FEM":False,
        "use_J":{"x": False, "y": False, "z":False}, 
        "use_rho":False, 
        "use_field": {"Ex":True, "Ey":False, "Ez":False, "Bx":False, "By":False, "Bz":False},
        "material":"vacuum", "gamma":0}
    ]
    k = open("data/FePt_bilayer-open/flux.out")
    spin_current = np.loadtxt("data/FePt_bilayer-open/flux.out") # 
    #print(spin_current.shape)
    sim = Simulation(spin_current, sim_params, med_params)
    sim.run()
    # print(sim.mediums[2].arrays["Ex"].shape)


