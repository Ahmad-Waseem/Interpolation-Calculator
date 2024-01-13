import numpy as np
import pandas as pd
import sympy as sp
from Newton_Backward import *
from Newton_Forward import *
from Gauss_Forward import *
from Gauss_Backward import *
from Bessel import *
from math import cos, log, sin, tan, exp, sqrt, ceil


import sympy as sp

# Usage example:
# df, delta_ys, poly, xnr, hr = newton_backward(df, num_of_entries, Xr)








num_of_entries = 10  #int(input("num of entries"))
x = [10,20,30,40,50,60,70,80,90,100]
y = [15.06223267,133.8123251,446.5988026,1024.771296,1911.199085,3113.597109,4600.245347,6298.387593,8095.50019,9843.47236]
#x = [2,3,4,5,6,7,8,9,10,11]
#y=[0.2621161459263087,1.2460441994649174,-0.20668693631404967,0.47580625413874394,0.07075429163087452,0.04638345169472506,0.006556352169799097,-0.004406808151794405,-0.0006363981162977492,0.0004101094814096555]
#variable to find
#x = [1,2,3,4,5,6,7,8,9,10]
#y = [-2.0001,-2.002,-0.012,3.96,9.904,17.802,27.63,39.37,53.002,68.480]
#

#x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

#y = [((exp(ex)* sin(ex)) / (tan(ex)*ex**3)) for ex in x]
#y = [1.263563,	1.991017,	1.521514,	0.693985,	0.504434,	0.8557,	1.108513,	0.956427,	0.578117,	0.270308]

Xr = 85
ex = Xr
org = ((exp(ex)* sin(ex)) / (tan(ex)*ex**3))
#org = 1.0685757
"""
for i in range(num_of_entries):
  val = float(input("enter the value of x"+ str(i)+ ": "))
  x.append(val)
  val = float(input("enter the value of y"+ str(i)+ ": "))
  y.append(val)
"""

data = {'x': x, 'delta_y_0': y}

# Create a Pandas DataFrame from the data dictionary
df = pd.DataFrame(data)

# Calculate and add the differences between consecutive rows
#
for zoro in range(3,4):
    indicator = zoro  #0 si backward, 1 is forward

    # Calculate and add the differences between consecutive rows in reverse order
    


    print(df)

    print(deltas)
    print(len(poly))
    # #print("\n\nFormula:")
    # polx = sp.latex(poly[-1])
    # #print(polx)
    # #print("\n\n Readable Formula:")
    # sp.pprint((poly[-1]))

    #print(sp.expand(poly))

