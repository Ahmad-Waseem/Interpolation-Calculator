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
        index = (df['x'] >= Xr).idxmax()               #To get id of P=0
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