import numpy as np

class Simulation:
    def __init__(self, spin_current, sim_params, medium_params):
        self.spin_current = spin_current
        self.sim_p = sim_params
        self.med_p = medium_params
        self.time = self.format_time()
        self.mediums = self.get_mediums()

    def format_time(self):
        dt = self.sim_p["dt"]
        time = np.arange(0, self.spin_current.shape[0] + 1, dt)
        # print(time)
        return time
    
    def get_mediums(self):
        N_time = len(self.time)
        mediums = []
        for medium in self.med_p:
            mediums.append(Medium(medium, self.sim_p, N_time))
        return mediums
    
    def step(self):
        dt = self.sim_p["dt"]
        
            
class Medium:
    def __init__(self, params, sim_params, N_time):
        self.params = params
        self.use_FEM = self.params['FEM']
        self.arrays = self.get_arrays(sim_params, N_time)
    
    def get_arrays(self, sim_p, N_time):
        dx, dy, dz = sim_p["dx"], sim_p["dy"], sim_p["dz"]
        N_points = lambda coord: int((self.params[coord][1] - self.params[coord][0])/sim_p[f"d{coord}"] + 1 - int(self.use_FEM))
        # use_FEM = self.params["FEM"]
        Nx = N_points("x")
        Ny = N_points("y")
        Nz = N_points("z")
        arrays = {}
        for coord in self.params["use_current"]:
            if self.params["use_current"][coord]:
                arrays[f"J{coord}"] = np.zeros((N_time, Nx,Ny,Nz))
            else:
                arrays[f"J{coord}"] = None
        if self.params["use_dist"]:
            arrays["rho"] = np.zeros((N_time, Nx,Ny,Nz))
        else:
            arrays["rho"] = None
        for coord in ["x", "y", "z"]:
            if self.params["use_field"]:
                arrays[f"E{coord}"] = np.zeros((N_time, Nx,Ny,Nz))
                arrays[f"B{coord}"] = np.zeros((N_time, Nx,Ny,Nz))
            else:
                arrays[f"E{coord}"] = None
                arrays[f"B{coord}"] = None
        
        return arrays

if __name__ == "__main__":
    jef_terms = {"rho": False, "rho_t": False, "J": False, "J_t": True} # Which terms from Jefimenko's equations to include
    exc_region = {"shape": "circle", "radius": 1000} # The region where the laser pulse hits the plates
    sim_params = {"dx":10, "dy":10, "dz":1, "dt": 1, # time is fs, length is nm
                "jef_terms": jef_terms, "exc_region": exc_region} 
    med_params = [
        {"type": "Mag",
        "x": (-2000, 2000), "y":(-2000, 2000), "z": (0,5), "FEM":True,
        "use_current":{"x": False, "y": False, "z":False}, 
        "use_dist":False, 
        "use_field": False,
        "material":"Fe", "gamma":0},

        {"type": "nonMag",
        "x": (-2000, 2000), "y":(-2000, 2000), "z": (5,8), "FEM":True,
        "use_current":{"x": True, "y": False, "z":False}, 
        "use_dist":False, 
        "use_field": False,
        "material":"Pt", "gamma":0.068},

        {"type":"vacuum",
        "x": (0, 0), "y":(0, 0), "z": (9, 10000), "FEM":False,
        "use_current":{"x": True, "y": False, "z":False}, 
        "use_dist":False, 
        "use_field": True,
        "material":"vacuum", "gamma":0}
    ]
    k = open("data/FePt_bilayer-open/flux.out")
    spin_current = np.loadtxt("data/FePt_bilayer-open/flux.out")
    print(spin_current.shape)
    sim = Simulation(spin_current, sim_params, med_params)


