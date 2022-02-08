"""This module contains functions that can help you manipulate many files
in your directories at once. It can check if folders contain a specific type
of file, and also delete files based on some critteria. This is useful for validating that 
all Fortran simulations went as expected, ie no simulation was corrupted"""

import os
import numpy as np
import shutil

def check_folders(folder_rootdir, filename):
    """This function checks all subdirectory of a directory if they contain a particular file

    If structure is:

    - folder_rootdir:
        - file1
        - simulation_results_1
            file11
            file12
            - some_folder11
                - file111
                - file112
        - ssimulation_results_2:
            - file21
    
    Then this function will check file11, file12 and file21

    The idea is that you should have a directory containing many folders with simulation results in each

    Arguments
    ---------
        folder_rootdir: str
            The root which we want to search throug
        
        filename: str
            The name of the file we look for. Could for example be 'flux.out'
    
    Returns
    -------
        corrupted_folders: List[str]
            List of directories where the file we looked for did not exist
        
        N_folders: int
            The number of corrupted folders (equivalent to len(corrupted_folders))

    """

    
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

def check_flux(folder_rootdir, expected_flux_time):
    """This function checks that the flux.out files in subdirectories of a root directory 
    are complete, i.e. that they contain all the desired timepoints. 
    
    In the Fortran simulations, sometimes you would get back a flux.out dile which was NOT complete.
    Tthis function lets you check all subdirectories od a directory to determine which that contain 
    complete/incomplete flux.out file
    
    This function can easily be copied and/or modified to look for other files

    The search sturcture is the same as in check_folders()
    
    Arguments
    ---------

        folder_rootdir: str
            The directory we want to search
        
        expected_flux_time: int
            The expected number of timesteps that the flux.out file should contain
    
    Returns
    -------
        complete_folders: List[str]
            List of the folder directories where the flux.out file looks correct
        
        incomplete_folders: List[str]
            List of the folder directories where the flux.out file looks incorrect
    
    """

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
    """
    USE THIS FUNCTION WITH CAUTION!!!

    This function recursively deletes all files of a specific type in a root directory

    Arguments
    ---------

        parent_dir: str
            The root directory from where we want to delete everything
        
        file_names: List[str]
            List of file names or file extensions that we wish to delete
    
    Returns
    -------
        None
    """
    file_death_note = []
    for dirpath, dirs, files in os.walk(parent_dir):
        for f in files:
            for ext in file_names:
                if (len(f) >= len(ext)) and (f[-len(ext):] == ext):
                        file_death_note.append("/".join([dirpath,f]))
            
    for f in file_death_note:
        os.remove(f)


def delete_dirtype(parent_dir, dir_names:list):
    """Same as delete_filetype, but instead you specify a folder name"""
    dir_death_note = []
    for dirpath, dirs, files in os.walk(parent_dir):
        for d in dirs:
            for ext in dir_names:
                if d[-len(ext):] == ext:
                    dir_death_note.append("/".join([dirpath,d]))
    for d in dir_death_note:
        shutil.rmtree(d)