import numpy as np
import pandas as pd
import sympy as sp
from Newton_Backward import *
from Newton_Forward import *
from Gauss_Forward import *
from Gauss_Backward import *
from Bessel import *
from math import cos, log, sin, tan, exp, sqrt, ceil

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
    if indicator == 0:  #Newton Backward
        print("\n=================================>NEWTON BACKWARD\n")
        for i in range(1, num_of_entries):
            df['delta_y_' + str(i)] = df['delta_y_' + str(i - 1)].diff()  
        #    for j in range(num_of_entries - 1, i - 1, -1):
         #       df.at[j, 'delta_y_'+str(i)] = table.at[j-1, 'delta_y_'+str(i - 1)] - table.at[j, 'delta_y_'+str(i - 1)]
        
        # 
        #delta_ys = get_sigma_ys(df)

        delta_ys = [sp.Symbol('delta_y_'+str(i)) for i in range(num_of_entries)]
        X = sp.symbols('x')
        xn = sp.symbols('xn')
        h = sp.symbols('h')
        
        
        poly = generate_NB_polynomial(num_of_entries, X, xn, h, delta_ys)
        xnr = df['x'].iloc[-1]
        hr = (df['x'].iloc[1]) - (df['x'].iloc[0])
        print("x = ", Xr)
        print("x0 = ", xnr)
        print("h = ", h)

        
        delta_ys = get_sigma_ys(df)
    #CONTINUE FROM SYMBOLING NB

    elif indicator == 1:  #Newton Forward
        print("\n=================================>NEWTON FORWARD\n")
        for i in range(1, num_of_entries):
            col_name = 'delta_y_'+str(i)
            df[col_name] = df['delta_y_'+str(i-1)].diff()

        delta_ys = [sp.Symbol('delta_y_'+str(i)) for i in range(num_of_entries)]
        X = sp.symbols('x')
        xn = sp.symbols('xn')
        h = sp.symbols('h')
        poly = generate_NF_polynomial(num_of_entries, X, xn, h, delta_ys)
        xnr = df['x'].iloc[0]
        hr = (df['x'].iloc[1]) - (df['x'].iloc[0])
        print("x = ", Xr)
        print("x0 = ", xnr)
        print("h = ", hr)
        delta_ys = get_delta_ys(df)

    elif indicator == 2:   #gauss forward
        print("\n=================================>GAUSS FORWARD\n")
        index = (df['x'] >= Xr).idxmax()                #To get id of x0

        for i in range(1, num_of_entries):                #generate delta_y_s in Dataframe
          df['delta_y_' + str(i)] = df['delta_y_' + str(i - 1)].diff()

        #put original values
        delta_ys = get_GF_delta_ys(num_of_entries,df , index)
        num_of_entries = len(delta_ys)                    #to avoid redundant processing
        deltas = delta_ys
        delta_ys = [sp.Symbol('delta_y_'+str(i)) for i in range(num_of_entries)]
        
        X = sp.symbols('x')
        xn = sp.symbols('xn')
        h = sp.symbols('h')
        poly = generate_GF_polynomial(num_of_entries, X, xn, h, delta_ys)
      
      #put original values
        delta_ys = deltas
        num_of_entries = len(delta_ys)                    #to avoid redundant processing
        xnr = df['x'][index]    #getting x0
        hr = (df['x'].iloc[1]) - (df['x'].iloc[0])
        print("x = ", Xr)
        print("x0 = ", xnr)
        print("h = ", hr)

    elif indicator == 3:   #gauss Backward
        print("\n=================================>GAUSS BACKWARD\n")
        index = (df['x'] >= Xr).idxmax()               #To get id
        print(index)
        for i in range(1, num_of_entries):                #generate delta_y_s in Dataframe
          df['delta_y_' + str(i)] = df['delta_y_' + str(i - 1)].diff()

        delta_ys = get_GB_delta_ys(num_of_entries,df , index)
        num_of_entries = len(delta_ys)                    #to avoid redundant processing
        deltas = delta_ys
        delta_ys = [sp.Symbol('delta_y_'+str(i)) for i in range(num_of_entries)]
        X = sp.symbols('x')
        xn = sp.symbols('xn')
        h = sp.symbols('h')
        poly = generate_GB_polynomial(num_of_entries, X, xn, h, delta_ys)

      #put original values
        delta_ys = deltas
        num_of_entries = len(delta_ys)                    #to avoid redundant processing
        print(num_of_entries)
        xnr = df['x'][index]    #getting x0
        hr = (df['x'].iloc[1]) - (df['x'].iloc[0])
        print("x = ", Xr)
        print("x0 = ", xnr)
        print("h = ", hr)


    elif indicator == 4:   #BESSEL
        print("\n=================================>BESSEL\n")
        index = (ceil(num_of_entries/2)-1)               #To get index
        print(index)
        for i in range(1, num_of_entries):                #generate delta_y_s in Dataframe
          df['delta_y_' + str(i)] = df['delta_y_' + str(i - 1)].diff()
      
        delta_ys = get_BSL_delta_ys(num_of_entries, df , index)
     
        deltas = delta_ys
        ylen = len(delta_ys)
        print(ylen,"<=")
        delta_ys = [sp.Symbol('delta_y_'+str(i)) for i in range(ylen)]
        X = sp.symbols('x')
        xn = sp.symbols('x0')
        h = sp.symbols('h')
        poly = generate_BSL_polynomial(num_of_entries, X, xn, h, delta_ys)

      #put original values
        delta_ys = deltas
        xnr = df['x'][index]    #getting x0
        hr = (df['x'].iloc[1]) - (df['x'].iloc[0])
        print("x = ", Xr)
        print("x0 = ", xnr)
        print("h = ", hr)



    print(df)

    print(deltas)
    print(len(poly))
    # #print("\n\nFormula:")
    # polx = sp.latex(poly[-1])
    # #print(polx)
    # #print("\n\n Readable Formula:")
    # sp.pprint((poly[-1]))

    #print(sp.expand(poly))

    for k in range(1, num_of_entries):
        print("\n          --------------Equation of Degree", k,  "----------------------------")
        poly[k] = poly[k].subs([('delta_y_'+str(i), delta_ys[i]) for i in range(num_of_entries)])
        poly[k] = poly[k].subs(xn, xnr)
        poly[k] = poly[k].subs(h, hr)
        
        # print("\n\nEquation:")
        # polx = sp.latex(poly[k])
        # print(polx)


        print("\n\nReadable non simplified:")
        sp.pprint((poly[k]))
        print("\n\nReadable Simplified-Equation:")
        #polz =sp.latex(sp.expand(poly))

        #print(polz)
        #print("\n\nReadable of above:")
        sp.pprint(sp.expand(poly[k]))


        poly[k] = poly[k].subs(X, Xr)
        print("\nf("+str(Xr)+") = ", poly[k])
        print("Original Val = ", org)
        print("difference = ", org- poly[k])
