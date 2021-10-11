import os
import re

def build_input_params_py(z, param_original_path, param_target_path, spin_flux_path):
    with open(param_original_path, 'r') as f:
        lines = f.read().split("\n")        
    lines[14] = f'\t"x": (-1000, 1000), "y":(-1000, 1000), "z": (0,{z[0]}),'
    lines[21] = f'\t"x": (-1000, 1000), "y":(-1000, 1000), "z": ({z[0]}, {sum(z)}),'
    lines[2] = f'spin_flux_path = "{spin_flux_path}"'
    input_params_string = "\n".join(lines)
    input_params_file = open(param_target_path,"w")
    input_params_file.write(input_params_string)
    input_params_file.close()

def build_simulation_setup(parent_dir, new_parent_dir, param_original_path):
    for folder in os.listdir(parent_dir):
        folder_dir = "/".join([parent_dir, folder])
        spin_flux_path = "/".join([folder_dir, "flux.out"])
        nums = re.findall(r"[0-9]+", folder)
        z = (int(nums[0]),int(nums[1]))
        param_target_path = "/".join([folder_dir, "input_params.py"])
        new_folder_dir = "/".join([new_parent_dir, folder])
        try:
            os.mkdir(new_folder_dir)
        except:
            pass
        build_input_params_py(z, param_original_path, param_target_path, spin_flux_path)