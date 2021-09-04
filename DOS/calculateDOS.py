import numpy as np


Pt = np.loadtxt("DOS/Pt/DOS.dat")
Fe = np.loadtxt("DOS/Fe/DOS.dat")

print(Pt.shape)
print(Fe.shape)

i = 0 
for row in Pt:
    E = row[0]
    if E < 3:
        print(row.shape)
        i += 1
        if i > 5:
            raise ValueError()