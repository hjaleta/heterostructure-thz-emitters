import os
import numpy as np
import shutil

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

def check_flux(folder_rootdir, expected_flux_time):#, remove, newdir):

    complete = 0
    complete_folders = []
    incomplete_folders = []
    
    for folder in os.listdir(folder_rootdir):
        folderpath = "/".join([folder_rootdir,folder])
        for f in os.listdir(folderpath):
            if f == "flux.out":
                try:
                    arr = np.loadtxt(folder_rootdir + "/" + folder + "/" + f)
                    t = arr.shape[0]
                    if t >= expected_flux_time:
                        complete_folders.append(folder)
                        complete += 1
                    else:
                        incomplete_folders.append(folder)
                except:
                    incomplete_folders.append(folder)
    return complete_folders, incomplete_folders

def delete_filetype(parent_dir, file_names:list):
    file_death_note = []
    for dirpath, dirs, files in os.walk(parent_dir):
        for f in files:
            for ext in file_names:
                if f[-len(ext):] == ext:
                    file_death_note.append("/".join([dirpath,f]))
    for f in file_death_note:
        os.remove(f)

def delete_dirtype(parent_dir, dir_names:list):
    dir_death_note = []
    for dirpath, dirs, files in os.walk(parent_dir):
        for d in dirs:
            for ext in dir_names:
                if d[-len(ext):] == ext:
                    dir_death_note.append("/".join([dirpath,d]))
    for d in dir_death_note:
        shutil.rmtree(d)