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



