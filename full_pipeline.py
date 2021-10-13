from E_fieldSimulation.simulation import Simulation
from freq_analysis.broadband import Signal
from pipeline_tools.check_data import check_folders, check_flux
from pipeline_tools.build_simulation_input import build_simulation_setup, build_param_json, decode_json_params
import os

spin_source_folder = "test_pipeline/Spin Currents"
sim_folder = "test_pipeline/Simulation Results"

# spin_source_folder = "Full Simulation/Spin Current Data/FuPt-closed"
# sim_folder = "Full Simulation/Simulation Results/FuPt-closed"
json_source_path = "Full Simulation/input_params.json"

check_source_folders = True
build_sim_params = True
run_simulation = True

sim_time = 600
signal_params = {"padfactor":9, "n_interpol":1, "dt":1, "bb_fraction": 10}

if check_source_folders:
    complete, non_complete = check_flux(spin_source_folder, sim_time)
    if len(non_complete)>0:
        print(f"Folders {', '.join(non_complete)} don't have a valid flux.out file ")
    else:
        print(f"All {len(complete)} folders have a valid flux.out file")

if build_sim_params:
    build_simulation_setup(spin_source_folder, sim_folder, json_source_path, num_layers=2)

if run_simulation:
    N_simulations = len(os.listdir(sim_folder))
    for folder in os.listdir(sim_folder):
        source_folderpath = "/".join( [spin_source_folder, folder])
        result_folderpath = "/".join( [sim_folder, folder])
        param_path = "/".join([source_folderpath, "input_params.json"])
        sim_params, medium_params, vacuum_params, spin_flux = decode_json_params(param_path, sim_time)
        sim = Simulation(spin_flux, sim_params, medium_params, vacuum_params, name = folder)
        print(medium_params)
        sim.run(print_percent = 1)
        for E in sim.vacuum.E_fields:
            # print(E.Ex)
            signal = Signal(E.Ex, signal_params, E.z, sim.name)
            filepath = "/".join([result_folderpath, f"Signal z = {E.z:.2E} nm "])
            fourier_plot_path = filepath + "spectra.png"
            BW_plot_path = filepath + "BW.png"
            transient_plot_path = filepath + "transient.png"
            json_path = filepath + "data.json"
            signal.plot_BW(BW_plot_path)
            signal.plot_signal(transient_plot_path)
            signal.plot_fourier(fourier_plot_path)
            signal.export_json(json_path)
        
