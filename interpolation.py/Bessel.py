#Functions for Gauss Backward

import pandas as pd
import numpy as np
import sympy as sp

def get_BSL_delta_ys(entries, df, index):
    delta_ys = []

    for i in range(entries):
        
        if (i % 2 == 0):
            delta_ys.append((df['delta_y_'+ str(i)][index]))
            index += 1
            
            if index + 1 > entries:    #if entries are odd, then last term will get only 1 of that
                return delta_ys
            delta_ys.append((df['delta_y_'+ str(i)][index]))
            
        else:
            delta_ys.append((df['delta_y_'+ str(i)][index]))
            
        

    return delta_ys



def generate_BSL_polynomial(num_terms, x, x0, h, delta_ys):

    index = 0
    polynomial = [((delta_ys[index] + delta_ys[index+1]) /2)]
    index += 2
    ylen = len(delta_ys)

    for i in range(1, num_terms):
        n = (-(i // 2))#-1
        if (i%2==0):#true
            
            term = ((delta_ys[index] + delta_ys[index+1]) /2)
            term /= sp.factorial(i)
            index += 2
            for j in range(i):#0->2=>  x-x0-1 x-x0
                term *= ((x - x0 ) / h) + n
                n += 1
            
            if(index+1 > ylen):
                
                polynomial.append(polynomial[-1] + term)
                return polynomial
          
            


                        
        else:
            term = (((x-x0)/h) - 1/2)*delta_ys[index]#i=3
            index += 1
            term /= sp.factorial(i)
            for j in range(1,i):    #1->3=> x-x0
                term *= ((x - x0 ) / h) + n
                n += 1
            
            if(index > ylen):
                print("mere kartoot")
                polynomial.append(polynomial[-1] + term)
                return polynomial
         
 
        polynomial.append(polynomial[-1] + term)

        
        
        
        
        

    return polynomial