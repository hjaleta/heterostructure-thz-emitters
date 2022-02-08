""" 
This module performs post-processing of the electric field signal
we find from the Simulation module
Intended unit for time is femtosecond, ie dt = 1 means dt = 1 ns
If using other unit, the axis labels of generated plots will be wrong
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq
from matplotlib.patches import Polygon
import json

class Signal():
    """
    This class takes a Electric field signal and calculates the broadband properties of it
    It has methods to generate plots of the transient signal, the fourier spectra or the broadband regime.

    Attributes
    ----------
        z: float
            The distance to the origin where the signal is being observed
        
        name: str (optional)
            The name of the Signal 
        
        p: dict
            This dictionary hold some parameters of the Signal. Values corresponding to the respective key is:

                'dt': int or float (default 1)
                    The time resolution of the input signal

                'padfactor': int (default 9)
                    How much do we pad the transient signal with zeros? 
                    (This is done to increase the resolution in Fourier Space)
                
                'n_interpol': int (default 1)
                    How many interpolation points do we add between each data point in the E-field transient signal
                    (This is done to increase the width in Fourier space)
                
                'bb_fraction': int or float (default 10)
                    This value defines how we calculate the broadband, see report for precise definition
                    If this value is 10, then broadband region is defined by 1/10 = 0.1 = 10 % of maximal amplitude
                    If it is 100, then this value instead is 1/100 = 1 %               

        t_org: np.array
            Original time vector
        
        E_org: np.array
            Original E-field vector
        
        t: np.array
            Interpolated and padded time vector
        
        E: np.array
            Interpolated and padded E-field vector
        
        freq: np.array
            frequency vector
        
        E_f: np.array of complex numbers
            Fourier transformed Electric field
        
        E_f_abs: np.array
            The modulus values of self.E_f
        
        df: float
            The frequency resolution
        
        bandwidth: float
            The calculated bandwidth of the fourier spectrum
        
        bb_freq: tuple (i0, i1, f0, f1)
            The borders of the broadband region 
            The left border f0 = self.freq[i0]
            The right border f1 = self.freq[i1]
    """
    def __init__(self, E, z, params = {},  name = "None"):
        self.z = z
        self.name = name
        self.p = self.unpack_params(params)
        self.check_params()
        self.t_org, self.E_org = self.get_original_signal(E)
        self.t, self.E = self.interpolate()
        self.pad()
        self.freq, self.E_f, self.E_f_abs, self.df = self.get_fourier()
        self.bandwidth, self.bb_freq = self.get_broadband()
        
    def unpack_params(self, params):
        default_params = {
                        "padfactor":9, 
                        "n_interpol":1, 
                        "dt":1, 
                        "bb_fraction": 10
                        }
        for key in params:
            default_params[key] = params[key]
        return default_params
        
    
    def check_params(self):
        if type(self.p["n_interpol"]) != int or self.p["n_interpol"] < 0:
            raise ValueError(f"params value for key 'n_interpol' must be positive integer, but given value is {self.p['n_interpol']}")
        if self.p["padfactor"] < 0:
            raise ValueError(f"Padding can't be negative, but given value is {self.p['padfactor']}")
        if self.p["dt"] < 0:
            raise ValueError(f"Timestep can't be negative, but given value is {self.p['dt']}")
        if type(self.p['bb_fraction']) != int:
            raise ValueError("bb_fraction must be a list of positive numbers")

    def get_original_signal(self, E):
        t = np.linspace(0, len(E)*self.p["dt"], len(E))
        return t, E

    def interpolate(self):
        t_pol = np.linspace(0, len(self.E_org)*self.p["dt"], len(self.E_org)*(self.p["n_interpol"] + 1))
        E_pol = np.interp(t_pol, self.t_org, self.E_org)
        self.p["dt"] /= (1+self.p["n_interpol"])
        return t_pol, E_pol

    def pad(self):
        E = self.E.copy()
        N = len(E)
        pad_points = self.p["padfactor"] * N 
        i_max, i_min = np.argmax(E), np.argmin(E)
        i_mid = (i_max + i_min)//2
        left_pad = 0
        right_pad = 0
        peak_to_left = i_mid < N/2

        if pad_points < abs(N-2*i_mid):
            if peak_to_left:
                left_pad += pad_points
            else:
                right_pad += pad_points
            pad_points = 0
        else:
            if peak_to_left:
                pad_points -= (N - 2*i_mid)
                left_pad += (N - 2*i_mid)
            else:
                pad_points -= (2*i_mid - N)
                right_pad += (2*i_mid - N)

        right_pad += pad_points//2
        left_pad += pad_points//2
        self.E = np.concatenate((np.zeros(left_pad), E, np.zeros(right_pad)))
        self.t = np.linspace(0, len(self.E)*self.p["dt"], len(self.E))

    def get_fourier(self):
        E_f = rfft(self.E)
        E_f_abs = np.abs(E_f)
        freq = rfftfreq(len(self.E), d = self.p["dt"])*1000
        df = freq[1] - freq[0]
        return freq, E_f, E_f_abs, df

    def get_broadband(self):
        E_f_max = np.amax(self.E_f_abs)
        bb_index  = np.argwhere(self.E_f_abs > E_f_max/self.p['bb_fraction']).flatten()
        i_min, i_max =  bb_index[[1,-1]]
        fmin, fmax = self.freq[[i_min, i_max]]
        delta_f = fmax - fmin
        bb_freq = (i_min, i_max, fmin, fmax)
        return delta_f, bb_freq
        
    def plot_fourier(self, figpath = "None"):
        fig, ax = plt.subplots()
        ax.plot(self.freq, self.E_f_abs)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Frequency [THz]")
        ax.set_ylabel("E-field [N/C]")
        ax.set_title(f"Frequency Spectra - {self.name}")

        if figpath != "None":
            if figpath[-4:] != ".png":
                figpath += ".png"
            fig.savefig(figpath)

    def plot_signal(self, figpath = "None"):
        fig, ax = plt.subplots()
        ax.plot(self.t_org, self.E_org)
        ax.set_xlabel("Time [fs]")
        ax.set_ylabel("E-field [N/C]")
        ax.set_title(f"Transient Signal - {self.name}")

        if figpath != "None":
            if figpath[-4:] != ".png":
                figpath += ".png"
            fig.savefig(figpath)
        
    def plot_BW(self, figpath = None):

        f_margin = 10
        i_min, i_max = self.bb_freq[0], self.bb_freq[1]
        plot_min = i_min - f_margin if i_min > f_margin else 0
        plot_max = i_max + f_margin if (i_max + f_margin) < len(self.freq) else len(self.freq)

        fill_f = self.freq[i_min:i_max+1].copy()
        fill_E = self.E_f_abs[i_min:i_max+1].copy()
        low_bound = 1/self.p['bb_fraction']
        verts = [(fill_f[0], low_bound), *zip(fill_f, fill_E), (fill_f[-1], low_bound)]
        poly = Polygon(verts, facecolor = "springgreen")
        
        fig, ax = plt.subplots()
        ax.add_patch(poly)
        ax.plot(self.freq[plot_min:plot_max], self.E_f_abs[plot_min:plot_max])
        ax.set_xlabel("Frequency [THz]")
        ax.set_ylabel("E-field [N/C]")
        ax.set_title(f"Frequency Spectra - {self.name}")
        text_x = (self.freq[i_min] * 2 + self.freq[i_max]) / 3
        text_y = np.amax(self.E_f_abs)/3
        ax.text(text_x, text_y, f"Bandwidth: {self.bandwidth:.2f} THz", ha = "center", fontsize = 12)
    
        if isinstance (figpath, str):
            if figpath[-4:] != ".png":
                figpath += ".png"
            fig.savefig(figpath)

    def export_json(self, path):
        signal = list(self.E_org.copy())
        fourier = list(self.E_f_abs)
        freq = list(self.freq)
        broadband = self.bandwidth
        name = self.name
    
        json_dict = {
            "name": name,
            "signal":signal,
            "fourier":fourier,
            "freq": freq,
            "broadband":broadband
        }
        with open(path, "w") as write_file:
            json.dump(json_dict, write_file)
        

if __name__ == "__main__":
    dt = 0.1 # fs
    E = np.sin(np.arange(0,100,dt))
    E = E.flatten()
    params = {"padfactor":9, "n_interpol":1, "dt":1, "bb_fraction": 10}

    S = Signal(E,params, 10**8)
    S.plot_signal()
    S.plot_BW("ExamplePlot.jpg")


