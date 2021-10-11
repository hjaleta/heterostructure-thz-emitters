import numpy as np

spin_flux_path = "Results - Spin Current/FuPt-closed/W6-5-L60/flux.out"
spin_flux = np.loadtxt(spin_flux_path)

# Which terms from Jefimenko's equations to include
jef_terms = {"rho": False, "rho_t": False, "J": False, "J_t": True} 
# The region where the laser pulse hits the plates. Shape can be "circle" or "rectangle"
# for circle, give radius in nm, for rectangle give side lengths x and y
exc_region = {"shape": "circle", "radius": 500, "x_side": 1000, "y_side": 1000} 
sim_params = {"dx":100, "dy":100, "dz":1, "dt": 1, # time is fs, length is nm
        "jef_terms": jef_terms, "exc_region": exc_region} 
medium_params = [
    {"material":"Fe",
    "x": (-1000, 1000), "y":(-1000, 1000), "z": (0,6),
    "use_Jx": False, "use_Jy": False, "use_Jz": False, 
    "use_rho":False, 
    "use_field": False,
    "theta":0, "mu":0, "n":1},

    {"material":"Pt",
    "x": (-1000, 1000), "y":(-1000, 1000), "z": (6,11),
    "use_Jx": True, "use_Jy": False, "use_Jz": False, 
    "use_rho":False, 
    "use_field": False,
    "theta":0.068, "mu": 0, "n": 70}
    ]

z_points = np.array([10**8])

vacuum_params =  {
    "label":"vacuum",
    "z": z_points,
    "use_J":{"x": False, "y": False, "z":False}, 
    "use_rho":False, 
    "use_field": {"Ex":True, "Ey":False, "Ez":False, "Bx":False, "By":False, "Bz":False},
    "material":"vacuum", "theta":0
    }

