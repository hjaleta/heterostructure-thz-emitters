import numpy as np
# Which terms from Jefimenko's equations to include
jef_terms = {"rho": False, "rho_t": False, "J": False, "J_t": True} 
# The region where the laser pulse hits the plates. Shape can be "circle" or "rectangle"
# for circle, give radius in nm, for rectangle give side lengths x and y
exc_region = {"shape": "circle", "radius": 500, "x_side": 1000, "y_side": 1000} 
sim_params = {"dx":500, "dy":500, "dz":1, "dt": 1, # time is fs, length is nm
        "jef_terms": jef_terms, "exc_region": exc_region} 
medium_params = [
    {"label": "Mag",
    "x": (-1000, 1000), "y":(-1000, 1000), "z": (0,6),
    "use_Jx": False, "use_Jy": False, "use_Jz": False, 
    "use_rho":False, 
    "use_field": False,
    "material":"Fe", "theta":0},

    {"label": "nonMag",
    "x": (-1000, 1000), "y":(-1000, 1000), "z": (6,11),
    "use_Jx": True, "use_Jy": False, "use_Jz": False, 
    "use_rho":False, 
    "use_field": False,
    "material":"Pt", "theta":0.068}
    ]

# z_points = np.array([1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]) * 100 # distances z for E(z,t) ranging between 0.1 um and 10 cm
z_points = np.array([1, 100, 10000, 10**6]) * 100 # distances z for E(z,t) ranging between 0.1 um and 10 cm


vacuum_params =  {
    "label":"vacuum",
    "z": z_points,
    "use_J":{"x": False, "y": False, "z":False}, 
    "use_rho":False, 
    "use_field": {"Ex":True, "Ey":False, "Ez":False, "Bx":False, "By":False, "Bz":False},
    "material":"vacuum", "theta":0
    }

spin_flux_path = "Results - Spin Current/FuPt-closed/W6-5-L60/flux.out"
spin_flux = np.loadtxt(spin_flux_path)