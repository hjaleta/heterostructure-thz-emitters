import os
import numpy as np
import matplotlib.pyplot as plt
import shutil
import re

def check_folders(folder_rootdir, filename):
    N_folders = 0
    corrupted_folders = []

    for folder in os.listdir(folder_rootdir):
        N_folders += 1
        folderpath = folder_rootdir+folder
        has_file = False
        for file in os.listdir(folderpath):
            if file == filename:
                has_file = True
        if not has_file:
            corrupted_folders.append(folderpath)
    
    return corrupted_folders, N_folders

def check_flux(folder_rootdir):#, remove, newdir):

    complete = 0
    complete_folders = []
    incomplete_folders = []
    
    for folder in os.listdir(folder_rootdir):
        folderpath = folder_rootdir+folder
        for file in os.listdir(folderpath):
            if file == "flux.out":
                try:
                    arr = np.loadtxt(folder_rootdir + "/" + folder + "/" + file)
                    t = arr.shape[0]
                    if t >= 600:
                        complete_folders.append(folder)
                        complete += 1
                    else:
                        incomplete_folders.append(folder)
                except:
                    incomplete_folders.append(folder)
    return complete_folders, incomplete_folders

def build_input_params_py(z, param_original_path, param_target_path, spin_flux_path):
    with open(param_original_path, 'r') as f:
        lines = f.read().split("\n")        
    lines[14] = f'"x": (-1000, 1000), "y":(-1000, 1000), "z": (0,{z[0]}),'
    lines[21] = f'"x": (-1000, 1000), "y":(-1000, 1000), "z": ({z[0]}, {sum(z)}),'
    lines[2] = f'spin_flux_path = "{spin_flux_path}"'
    input_params_string = "\n".join(lines)
    input_params_file = open(param_target_path,"w")
    input_params_file.write(input_params_string)
    input_params_file.close()

def delete_filetype(parent_dir, file_extensions:list):
    all_filepaths = []
    i = 0
    j = 0
    for dirpath, dirs, files in os.walk(parent_dir):
        for file in files:
            for ext in file_extensions:
                if file[-len(ext):] == ext:
                    os.remove("/".join([dirpath,file]))

def build_simulation_setup(parent_dir, new_parent_dir, param_original_path):
    for folder in os.listdir(parent_dir):
        folder_dir = "/".join([parent_dir, folder])
        spin_flux_path = "/".join([folder_dir, "flux.out"])
        nums = re.findall(r"[0-9]+", folderpath)
        z = (int(nums[0]),int(nums[1]))
        param_target_path = "/".join([folder_dir, "input_params.py"])
        new_folder_dir = "/".join([new_parent_dir, folder])
        os.mkdir(new_folder_dir)
        build_input_params_py(z, param_original_path, param_target_path, spin_flux_path)


folder_dir = "Results - Spin current/FuPt-closed"
new_folder_dir = "Simulation Input/FuPt-closed/"
param_original_path = "E_fieldSimulation/input_params.py"

build_input_params_py((2,4), param_original_path, "oneparamtest.py", "Results - Spin current/FuPt-closed/W2-5-L60/flux.out")


# In[ ]:




