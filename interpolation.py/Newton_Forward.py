#Functions for Newton Forward
import sympy as sp
import pandas as pd
import numpy as np


def get_delta_ys(df):
  # Initialize a list for delta_y values
  delta_ys = [df['delta_y_0'][0]]

  # Loop through delta_y columns and fetch the desired terms
  for i in range(1, len(df.columns) - 1):  # Exclude 'x' and 'y' columns
    column_name = 'delta_y_' + str(i)
    delta_y_value = df[column_name][i]
    delta_ys.append(delta_y_value)

  return delta_ys


#-------------------------------------------------------


def generate_NF_polynomial(entries, x, xn, h, delta_ys):
    # Initialize the polynomial
    polynomial = [delta_ys[0]]

    # Calculate and add the subsequent terms to the polynomial
    for i in range(1, entries):
        
        term = delta_ys[i] / sp.factorial(i)
        for j in range(i):
            term *= ((x - xn) / h) -j
        polynomial.append(polynomial[-1] + term)
    return polynomial


