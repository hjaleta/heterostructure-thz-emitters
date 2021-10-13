import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq
from matplotlib.patches import Polygon
import pandas as pd
import json

class Signal():
    def __init__(self, E, params, z, name = "None"):
        self.z = z
        self.name = name
        self.unpack_params(params)
        self.check_params()
        self.t_org, self.E_org = self.get_original_signal(E)
        self.t, self.E = self.interpolate()
        self.pad()
        self.freq, self.E_f, self.E_f_abs, self.df = self.get_fourier()
        self.bandwidth, self.bb_freq = self.get_broadband()
        
    def unpack_params(self, params):
        for key in params:
            exec(f"self.{key} = params[key]")
    
    def check_params(self):
        if type(self.n_interpol) != int:
            raise ValueError(f"Number of interpolation points must be type int, but given type is {type(self.n_interpol)}")
        elif self.n_interpol < 0:
             raise ValueError("Number of interpolation points cant be negative")
        if self.padfactor < 0:
            raise ValueError(f"Padding can't be negative, but given value is {self.P}")
        if self.dt < 0:
            raise ValueError(f"Timestep can't be negative, but given value is {self.dt}")
        if type(self.bb_fraction) != int:
            raise ValueError("bb_fraction must be a list of positive numbers")

    def get_original_signal(self, E):
        t = np.linspace(0, len(E)*self.dt, len(E))
        return t, E

    def interpolate(self):
        t_pol = np.linspace(0, len(self.E_org)*self.dt, len(self.E_org)*(self.n_interpol +1))
        E_pol = np.interp(t_pol, self.t_org, self.E_org)
        self.dt /= (1+self.n_interpol)
        return t_pol, E_pol

    def pad(self):
        E = self.E.copy()
        N = len(E)
        pad_points = self.padfactor * N
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
        self.t = np.linspace(0, len(self.E)*self.dt, len(self.E))

    def get_fourier(self):
        E_f = rfft(self.E)
        E_f_abs = np.abs(E_f)
        freq = rfftfreq(len(self.E), d = self.dt)*1000
        df = freq[1] - freq[0]
        return freq, E_f, E_f_abs, df

    def get_broadband(self):
        E_f_max = np.amax(self.E_f_abs)
        bb_index  = np.argwhere(self.E_f_abs > E_f_max/self.bb_fraction).flatten()
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
        ax.set_title("Transient Signal - {self.name}")

        if figpath != "None":
            if figpath[-4:] != ".png":
                figpath += ".png"
            fig.savefig(figpath)
        
    def plot_BW(self, figpath = "None"):

        f_margin = 10
        i_min, i_max = self.bb_freq[0], self.bb_freq[1]
        plot_min = i_min - f_margin if i_min > f_margin else 0
        plot_max = i_max + f_margin if (i_max + f_margin) < len(self.freq) else len(self.freq)

        fill_f = self.freq[i_min:i_max+1].copy()
        fill_E = self.E_f_abs[i_min:i_max+1].copy()
        low_bound = 1/self.bb_fraction
        verts = [(fill_f[0], low_bound), *zip(fill_f, fill_E), (fill_f[-1], low_bound)]
        poly = Polygon(verts, facecolor = "springgreen")
        
        fig, ax = plt.subplots()
        ax.add_patch(poly)
        ax.plot(self.freq[plot_min:plot_max], self.E_f_abs[plot_min:plot_max])
        ax.set_xlabel("Frequency [THz]")
        ax.set_ylabel("E-field [N/C]")
        ax.set_title("Frequency Spectra - {self.name}")
        text_x = (self.freq[i_min] * 2 + self.freq[i_max]) / 3
        text_y = np.amax(self.E_f_abs)/3
        ax.text(text_x, text_y, f"Bandwidth: {self.bandwidth:.2f} THz", ha = "center", fontsize = 12)
    
        if figpath != "None":
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
    E = np.loadtxt("../testsave copy.out")
    E = E[:,-1].flatten()
    params = {"padfactor":9, "n_interpol":1, "dt":1, "bb_fraction": 10}

    S = Signal(E,params)
    # S.plot_fourier()
    S.plot_signal()
    S.plot_BW()
    print(S.freq[:5], S.freq[-5:])



