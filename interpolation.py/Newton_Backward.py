#Functions for Newton Backward

import pandas as pd
import numpy as np
import sympy as sp

def get_sigma_ys(df):
  # Initialize a list for last non-null delta_y values
  delta_ys = [df['delta_y_0'].iloc[-1]]  # Initialize with the last 'y' value

  p = len(df.columns) - 2  # Index for traversing 'delta_y' columns in reverse order

  # Loop through delta_y columns in reverse and fetch the last non-null values
  for i in range(1, len(df.columns) - 1):
    column_name = 'delta_y_' + str(i)
    delta_y_value = df[column_name].dropna().iloc[-1]
    delta_ys.append(delta_y_value)
    p -= 1

  return delta_ys



#-------------------------------------------------------


def generate_NB_polynomial(entries, x, xn, h, delta_ys):
    # Initialize the polynomial
    polynomial = [delta_ys[0]]

    # Calculate and add the subsequent terms to the polynomial
    for i in range(1, entries):
        term = delta_ys[i] / sp.factorial(i)
        for j in range(i):
            term *= ((x - xn) / h) + j
        polynomial.append(polynomial[-1] + term)
    return polynomial


#----------------------------- Calculations

def newton_backward(df, num_of_entries, Xr):
    print("\n=================================>NEWTON BACKWARD\n")
    for i in range(1, num_of_entries):
        df['∇y_' + str(i)] = df['∇y_' + str(i - 1)].diff()  

    delta_ys = [sp.Symbol('∆y_'+str(i)) for i in range(num_of_entries)]  # Changed to triangular delta symbol (∆)
    X = sp.symbols('X')
    xn = sp.symbols('Xn')
    h = sp.symbols('h')

    poly = generate_NB_polynomial(num_of_entries, X, xn, h, delta_ys)

    xnr = df['x'].iloc[-1]
    hr = (df['x'].iloc[1]) - (df['x'].iloc[0])
    print("X = ", Xr)
    print("Xn = ", xnr)
    print("h = ", h)

    delta_ys = get_sigma_ys(df) 

    return df, delta_ys, poly, xnr, hr
#Usage: newton_backward(pd, np-ndarray, X to find-double)



#------------------- Printing
def print_GB(poly, delta_ys, xnr, hr, Xr, org, num_of_entries, degree):
    for k in range(1, num_of_entries):
        print("\n          --------------Equation of Degree", k,  "----------------------------")
        poly[k] = poly[k].subs([(f'∇y_{i}', delta_ys[i]) for i in range(num_of_entries)])
        poly[k] = poly[k].subs('xn', xnr)
        poly[k] = poly[k].subs('h', hr)

        print("\n\nReadable non-simplified:")
        sp.pprint((poly[k]))
        print("\n\nReadable Simplified-Equation:")
        sp.pprint(sp.expand(poly[k]))

        poly[k] = poly[k].subs('x', Xr)
        print("\nf(", Xr, ") =", poly[k])
        print("Original Val =", org)
        print("Difference =", org - poly[k])

# Usage example: print_GB(poly, delta_ys, xnr, hr, Xr, org, num_of_entries)
