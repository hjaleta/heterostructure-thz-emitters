"""
This model computes the best linear regression fit to the results
The purpose is to investigate the correlation between the features (widths and laser pulse)

"""

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing

def fitResults(df, include_crossterms=False):
    arr = df.to_numpy(dtype=np.float64)
    Y = arr[:,0].reshape(-1,1)
    X = arr[:,1:]
    if include_crossterms:
        cross_terms = []
        cross_terms.append(np.array(X[:,0]*X[:,1]))
        cross_terms.append(np.array(X[:,0]*X[:,2]))
        cross_terms.append(np.array(X[:,1]*X[:,2]))
        ct = np.array(cross_terms).T
        X = np.concatenate((X,ct), axis=1)

    scaler = preprocessing.StandardScaler().fit(X)
    X_scaled = scaler.transform(X)
    linear_model = LinearRegression().fit(X_scaled, Y)
    coeff_dict = {feature:coef for feature, coef in zip(list(df.columns)[1:], list(linear_model.coef_)[0]) }
    coeff_dict["Intercept"] = float(linear_model.intercept_)
    return coeff_dict

if __name__ == "__main__":
    from AnalyseResults.ResultToDF import getDF
    result_dir = "Full Simulation/Simulation Results/FuPt-closed/"
    df = getDF(result_dir)
    coeff_dict = fitResults(df, False)
    print(coeff_dict)