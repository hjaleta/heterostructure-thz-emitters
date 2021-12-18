import json
import os
import sys
cwd = os.getcwd()
sys.path.append(cwd[:-15])
import pandas as pd
import re

def getDF(result_parent_dir, sort_df = False):
    cwd = os.getcwd()
    all_results_path = cwd + "/" + result_parent_dir

    data_dict = {"Signal Bandwidth [THz]": [], "Width Iron [nm]":[], "Width Platinum [nm]":[], "Laser Pulse Time [fs]": []}

    num_pattern = r"W(\d)-(\d)-L(\d*)"
    print(all_results_path)
    for dirpath, dirs, files in os.walk(all_results_path):
        for f in files:
            if f[-4:] == "json":
                f_path = dirpath + "/" + f
                
                with open(f_path, "r") as j_object:
                    results = json.load(j_object)
                m = re.search(num_pattern, results["name"])
                W1, W2, L = m.group(1,2,3)
                BW = results["broadband"]
                data_dict["Signal Bandwidth [THz]"].append(BW)
                data_dict["Width Iron [nm]"].append(W1)
                data_dict["Width Platinum [nm]"].append(W2)
                data_dict["Laser Pulse Time [fs]"].append(L)
                
    df = pd.DataFrame.from_dict(data_dict)
    if sort_df:
        df = df.sort_values(["Signal Bandwidth [THz]"], ignore_index=True, ascending=False)
        df.index = range(1, len(df)+1)
    return df

def export_df(df, path):
    df.to_csv(path, float_format= '%.1f')

if __name__ == "__main__":
    # result_parent_dir = "Full Simulation/Simulation Results/FuPt-closed/" 
    # export_dir = "/home/hjaleta/spin_project/AnalyseResults/Sorted Table of Bandwidths2.csv"
    # df = getDF(result_parent_dir,True)
    # export_df(df, export_dir)
    pass