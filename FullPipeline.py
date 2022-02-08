""" 
This script runs several E-Field simulations. It obviously requires some folders which hold results from the
Fortran simulations
"""

from EFieldSimulation.Simulation import Simulation
from EFieldSimulation.Frequency.Signal import Signal
from PipelineTools.FileManagement import check_folders, check_flux
from PipelineTools.build_simulation_input import build_simulation_setup, build_param_json, decode_json_params
import os

# Folder with subfolders with results from Fortran

# spin_source_folder = "FortranSimulation/Fortran Results/FuPt-closed"
spin_source_folder = "FortranSimulation/Fortran Results/test_results"

# Folder to save the subfolders with simulation results
sim_folder = "EFieldSimulation/Simulation Results_2/FuPt-closed"

# Template file for input parameters 
json_source_path = "EFieldSimulation/input_params.json"

# Check that Fortran result folders are not corrupted
check_source_folders = True

# Construct json-files with simulation parameters
build_sim_params = True

# Run all simulations
run_simulation = True

# Number of timesteps
sim_time = 600

# Parameters for signal post processing
signal_params = {"padfactor":9, "n_interpol":1, "dt":1, "bb_fraction": 10}

def main():
    if check_source_folders: # Check validity of flux.out files
        complete, non_complete = check_flux(spin_source_folder, sim_time)
        if len(non_complete)>0:
            print(f"Folders {', '.join(non_complete)} don't have a valid flux.out file ")
        else:
            print(f"All {len(complete)} folders have a valid flux.out file")

    if build_sim_params: # Construct the necessary json files
        build_simulation_setup(spin_source_folder, sim_folder, json_source_path, num_layers=2)

    if run_simulation: # Run simulation
        N_simulations = len(os.listdir(sim_folder)) # Number of simulations
        for f_i, folder in enumerate(os.listdir(sim_folder)): # Iterate iver result folders

            try: # Try to perform the simulation
                # Define necessary paths
                source_folderpath = "/".join( [spin_source_folder, folder])
                result_folderpath = "/".join( [sim_folder, folder])
                param_path = "/".join([source_folderpath, "input_params.json"])

                # Extract simulation parameters from json file
                sim_params, medium_params, vacuum_params, spin_flux = decode_json_params(param_path, sim_time)

                # Create simulation object
                sim = Simulation(spin_flux, sim_params, medium_params, vacuum_params, name = folder)
                # print(medium_params)
                # print(signal_params)
                # Run simulation
                sim.run()

                for E in sim.vacuum.E_fields: # Iterate over output signals

                    # Create Signal object
                    print("signal length ", len(E.Ex))
                    signal = Signal(E.Ex,  E.z, signal_params, sim.name)

                    # Generate filepaths for results
                    filepath = "/".join([result_folderpath, f"Signal z = {E.z:.2E} nm "])
                    fourier_plot_path = filepath + "spectra.png"
                    BW_plot_path = filepath + "BW.png"
                    transient_plot_path = filepath + "transient.png"
                    json_path = filepath + "data.json"

                    # Save plots and data
                    signal.plot_BW(BW_plot_path)
                    signal.plot_signal(transient_plot_path)
                    signal.plot_fourier(fourier_plot_path)
                    signal.export_json(json_path)

                print(f"Simulation {f_i + 1} of {N_simulations} completed")

            except: # Declare failure
                print(f"Simulation {f_i + 1} of {N_simulations} did not complete")

if __name__ == "__main__":
    main()