from E_fieldSimulation.simulation import Simulation
from freq_analysis.broadband import Signal
from pipeline_tools.check_data import check_folders, check_flux
from pipeline_tools.build_simulation_input import build_input_params_py, build_simulation_setup
import os
import sys

spin_source_folder = "test_pipeline/Spin Currents"
sim_folder = "test_pipeline/Simulation Results"

check_source_folders = True
build_sim_params = True
run_simulation = False

sim_time = 600
signal_params = {"padfactor":9, "n_interpol":1, "dt":1, "bb_fraction": 10}

if check_source_folders:
    complete, non_complete = check_flux(spin_source_folder, sim_time)
    if len(non_complete)>0:
        print(f"Folders {', '.join(non_complete)} don't have a valid flux.out file ")
    else:
        print(f"All {len(complete)} folders have a valid flux.out file")

if build_sim_params:
    pass

if run_simulation:
    
    N_simulations = len(os.listdir(sim_folder))
    for folder in os.listdir(sim_folder):
        folderpath = "/".join( [sim_folder, folder])
        sys.path.append(folderpath)
        from input_params import sim_params, medium_params, vacuum_params, spin_flux
        sys.path.remove(folderpath)
        sim = Simulation(spin_flux, sim_params, medium_params, vacuum_params, name = folder)
        sim.run(print_progress=True)
        for E in sim.vacuum.E_fields:
            signal = Signal(E.Ex, signal_params)
            BW_plot_path = "/".join([folderpath, f"Signal z = {E.z:.2E} BW.png" ])