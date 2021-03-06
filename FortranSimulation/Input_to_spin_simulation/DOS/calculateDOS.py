""" This script contains functions that calculate the electron / (cell*spin*energy evel)
    density, a sort of density of states 

    The input data must be compatible with numpy. It uses txt- or dat-files but could easily
    be modified to open other file types. The input data should have shape (4, N) where N
    is the number of rows. The four columns should be [E, N(up), N(down), N(tot)]. 
    E is the energy in eV, N(up) is the number of spin up electrons per energy level and unit cell.
    N(down) is the same as N(up) but for spin up, N(tot) is N(up) + N(down)
"""

import numpy as np
     
def trapz(f, x):
    # Trapezoidal integration
    # f is vector of functional values f(x)
    vals = []
    for a, b, x0, x1 in zip(f[:-1], f[1:], x[:-1], x[1:]):
        dx = x1 - x0
        val = (a+b)*dx/2
        vals.append(val)
    vals_arr = np.array(vals)
    return vals_arr

def calculateDOS(data_path, filename, E_0, N_channels, dE, scaling_factor = 1):
    """ 
    Integrates the density of states from the data files

    Args:
        data_path (string): path to the data file
        filename (string): path of the output file 
        E_0 (float): Lower limit of the first energy channel to be calculatetd
        N_channels (int): Number of channels to be calculated
        dE (float): width of energy channels
        scaling_factor (float): Factor to rescale result. For example from unit cell to nm^2. Defaults to 1
    """

    try: # Load data 
        data = np.loadtxt(data_path)
    except:
        raise ValueError("Provided DOS data not compatible")
    if data.shape[1] != 4:
        raise ValueError("Provided DOS data not compatible")
    if type(N_channels) != int or N_channels < 1:
        raise ValueError("N_channels must be positive integer")

    E_disc = np.arange(E_0, E_0 + (N_channels+1)*dE, dE)
    spin_up_pdf = np.interp(E_disc, data[:,0].flatten(), data[:,1].flatten())
    spin_down_pdf = np.interp(E_disc, data[:,0].flatten(), data[:,2].flatten())
    spin_up = trapz(spin_up_pdf, E_disc).reshape((1,-1))
    spin_down = trapz(spin_down_pdf, E_disc).reshape((1,-1))
    spin_tot = spin_up + spin_down
    all_spin = np.concatenate((spin_up, spin_down, spin_tot))*scaling_factor
    all_spin = all_spin.T
    header_str = f"Number of electrons per channel\nSpin up\tSpin down\tTotal"
    np.savetxt(filename, all_spin, fmt='%1.3f', header = header_str)


if __name__ == "__main__":

    Pt_datapath = "/home/hjaleta/spin_project/DOS/Pt/DOS.dat"
    Pt_filepath = "/home/hjaleta/spin_project/DOS/Pt/Energy channels Pt.txt"
    Pt_unit_cell_volume_inv = 1000/60.43

    Fe_datapath = "/home/hjaleta/spin_project/DOS/Fe/DOS.dat"
    Fe_filepath = "/home/hjaleta/spin_project/DOS/Fe/Energy channels Fe.txt"
    Fe_unit_cell_volume_inv = 1000/23.55

    E_0 = 0
    N_channels = 12
    dE = 0.125

    calculateDOS(Pt_datapath, Pt_filepath, E_0, N_channels, dE, Pt_unit_cell_volume_inv)
    calculateDOS(Fe_datapath, Fe_filepath, E_0, N_channels, dE, Fe_unit_cell_volume_inv)