"""
This model computes the best linear regression fit to the results
The purpose is to investigate the correlation between the features (widths and laser pulse)
and the output (broadband)
"""

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
import os
import json
import sys
import pandas as pd
import re


def getDF(result_parent_dir, sort_df = False):
    """
    Recursively takes all files from a directory (top to bottom), finds all 
    json result files and extract the 2 widths, the length of the laser pulse and the broadband
    from each file
    Args:
        result_parent_dir: str
            The directory where we save the results
        
        sort_df: bool
            Determines whether the df should be sorted after broadband

    Returns:
        df: pandas.DataFrame
            The DataFrame with the results
        
    """

    # Get the location of the json files
    cwd = os.getcwd()
    all_results_path = cwd + "/" + result_parent_dir

    # Initiate dict with data
    data_dict = {"Signal Bandwidth [THz]": [], "Width Iron [nm]":[], "Width Platinum [nm]":[], "Laser Pulse Time [fs]": []}

    # Regex patern to find features
    num_pattern = r"W(\d)-(\d)-L(\d*)"
    for dirpath, dirs, files in os.walk(all_results_path): # Iterate over all subdirectories
        for f in files:
            if f[-4:] == "json":
                f_path = dirpath + "/" + f
                
                with open(f_path, "r") as j_object:
                    results = json.load(j_object)
                
                # Use regex pattern to find W1 W2 and L
                m = re.search(num_pattern, results["name"])
                W1, W2, L = m.group(1,2,3)

                # Get bandwidth
                BW = results["broadband"]

                # Put extracted dat a into dictionary
                data_dict["Signal Bandwidth [THz]"].append(BW)
                data_dict["Width Iron [nm]"].append(W1)
                data_dict["Width Platinum [nm]"].append(W2)
                data_dict["Laser Pulse Time [fs]"].append(L)
    
    # Create DataFrame from dictionary
    df = pd.DataFrame.from_dict(data_dict)

    if sort_df: # Sort dataframe after bandwidth
        df = df.sort_values(["Signal Bandwidth [THz]"], ignore_index=True, ascending=False)
        df.index = range(1, len(df)+1)
    return df

def export_df(df, path):
    """Save the DataFrame at the speciifed path"""
    df.to_csv(path, float_format= '%.1f')

def fitResults(df, include_crossterms=False):

    # Initialize data point
    arr = df.to_numpy(dtype=np.float64)
    Y = arr[:,0].reshape(-1,1)
    X = arr[:,1:]

    # If we want crosstemrms, initialise their data
    if include_crossterms:
        cross_terms = []
        cross_terms.append(np.array(X[:,0]*X[:,1]))
        cross_terms.append(np.array(X[:,0]*X[:,2]))
        cross_terms.append(np.array(X[:,1]*X[:,2]))
        ct = np.array(cross_terms).T
        X = np.concatenate((X,ct), axis=1)

    # Normalise data
    scaler = preprocessing.StandardScaler().fit(X)
    X_scaled = scaler.transform(X)

    # Fit data
    linear_model = LinearRegression().fit(X_scaled, Y)

    # Save cofficients in dictionary
    coeff_dict = {feature:coef for feature, coef in zip(list(df.columns)[1:], list(linear_model.coef_)[0]) }
    coeff_dict["Intercept"] = float(linear_model.intercept_)
    return coeff_dict

if __name__ == "__main__":
    from AnalyseResults.ResultToDF import getDF
    result_dir = "Full Simulation/Simulation Results/FuPt-closed/"
    df = getDF(result_dir)
    coeff_dict = fitResults(df, False)
    print(coeff_dict)