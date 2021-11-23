import os
import re
import json
import numpy as np

def build_param_json(source_path, target_path, z, spin_flux_path):
    with open(source_path, "r") as source_params:
        params = json.load(source_params)
    if len(z) != len(params["medium_params"]):
        raise ValueError("Faulty Z or json input")
    if target_path[-5:] != ".json":
        raise ValueError("target_path must be json")
    for p_i, pair in enumerate(z):
        params["medium_params"][p_i]["z"] = z[p_i]
    params["spin_flux_path"] = spin_flux_path
    j_object = json.dumps(params, indent=4)
    with open (target_path, "w") as json_param_file:
        json_param_file.write(j_object)

def decode_json_params(json_param_path, t_max):
    with open(json_param_path, "r") as json_params:
        params = json.load(json_params)
    spin_flux = np.loadtxt(params["spin_flux_path"])
    spin_flux = spin_flux[:t_max,:]
    sim_params = params["sim_params"]
    medium_params = params["medium_params"]
    vacuum_params = params["vacuum_params"]
    vacuum_params["z"] = np.array(vacuum_params["z"])
    return sim_params, medium_params, vacuum_params, spin_flux

def build_simulation_setup(parent_dir, new_parent_dir, json_source_path, num_layers=2):
    for folder in os.listdir(parent_dir):
        folder_dir = "/".join([parent_dir, folder])
        spin_flux_path = "/".join([folder_dir, "flux.out"])
        nums = re.findall(r"[0-9]+", folder)
        z_thick = [int(nums[i]) for i in range(num_layers)]
        z_min = 0
        medium_z_range = []
        for t in z_thick:
            medium_z_range.append((z_min, z_min + t))
            z_min = t
        json_target_path = "/".join([folder_dir, "input_params.json"])
        new_folder_dir = "/".join([new_parent_dir, folder])
        try:
            os.mkdir(new_folder_dir)
        except:
            pass
        build_param_json(json_source_path, json_target_path, medium_z_range, spin_flux_path)