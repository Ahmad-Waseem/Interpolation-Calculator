#Functions for Newton Backward

import pandas as pd
import numpy as np
import sympy as sp

def get_GF_delta_ys(num, df, index):
    entries = num
    inc = 0
    delta_ys = []

    for i in range(entries):
        
        x = df['delta_y_'+ str(i)][index + inc]
        if np.isnan(x):
            return delta_ys  # Return the list if 'x' is NaN
        delta_ys.append(x)
        if i % 2 == 0:
            inc += 1
        if index + inc >= entries:
            return delta_ys
        print("pass ", i)

    return delta_ys


def generate_GF_polynomial(num_terms, x, x0, h, delta_ys):

    polynomial = [delta_ys[0]]

    for i in range(1, num_terms):
        n = int(-(i / 2))
        term = delta_ys[i]/sp.factorial(i)
        for j in range(i):
            term *= ((x - x0 ) / h) + n
            n += 1
        
        polynomial.append(polynomial[-1] + term)

    return polynomial