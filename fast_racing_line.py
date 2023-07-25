

from casadi import Opti
from casadi import sumsqr
from casadi import diff
from casadi import mtimes

import numpy as np
from scipy.linalg import expm
from scipy.signal import cont2discrete
from matplotlib import pyplot as plt

from vehicle import *
from boostrtrees import RTree

from utility import *


def getAffineModel(Ux, K, ts, veh, afHat, arHat, counter):

    a = veh['a']; 
    b = veh['b']; 
    m = veh['m']; 
    Cf = veh['Cf']; 
    Cr = veh['Cr']; 
    g = veh['g']; 
    L = a+b; 
    Iz = veh['Iz'];
    FzF = veh['FzF']; 
    FzR = veh['FzR']; 
    muF = veh['muF']; 
    muR = veh['muR']; 
    xLA = veh['xLA']; 
    kLK = veh['kLK'];

    #k1 = kLK; k2 = xLA*kLK; #Note: In optimization, we don't want feedback acting.
    k1 = 0; 
    k2 = 0;

    FyFhat = tireforces(Cf, muF, muF, afHat, FzF, counter);
    FyRhat = tireforces(Cr, muR, muR, arHat, FzR, counter);  
    Cf_b = getLocalStiffness(afHat, Cf, muF, muF, FzF, counter);
    Cr_b = getLocalStiffness(arHat, Cr, muR, muR, FzR, counter);

    
    Ac = np.array([
                    [0, Ux, 0, Ux, 0],
                    [0, 0, 1, 0, 0], 
                    [-a*k1*Cf_b/Iz,   -a*k2*Cf_b/Iz,      (-a**2*Cf_b - b**2*Cr_b)/(Ux*Iz),   (b*Cr_b - a*Cf_b)/Iz, 0],
                    [-k1*Cf_b/(m*Ux),  -k2*Cf_b/(m*Ux),   (b*Cr_b - a*Cf_b)/(m*Ux**2)-1,     -(Cf_b+Cr_b)/(m*Ux), 0],
                    [0, 0, 1, 0, 0]
                  ]);

    Bc = np.array([[0], [0], [a*Cf/Iz], [Cf/(m*Ux)], [0]]);
     
    dc = np.array([[0], [-K*Ux], [(a*Cf*afHat - b*Cr*arHat)/Iz + (a*FyFhat-b*FyRhat)/Iz], [(Cf*afHat + Cr*arHat)/(m*Ux) + (FyFhat + FyRhat)/(m*Ux)], [0]]);
    
    Cc = np.ones((1,5))
    Dc = np.zeros((1,1))

    disc = cont2discrete((Ac,Bc, Cc, Dc), ts)

    return disc[0], disc[1], dc, Cf, Cr
    #return Ac, Bc, dc, Cf, Cr

def getAllSys(veh, Ux, K, ts):

    N = len(Ux);
    A = np.zeros((N, 5, 5));
    B = np.zeros((N,5,1));
    D = np.zeros((N,5,1));
    
    deltaFFW = [] 
    betaFFW = [] 
    
    aF = veh['alphaFrontTable']
    aR = veh['alphaRearTable']
    perc = 0
    
    print("-Computing system dynamic...")
    for i in range(0,N):
        
        perc = print_perc(perc, N, i)
    
        aFHat = force2alpha(veh['FyFtable'], veh['alphaFrontTable'], veh['b']/veh['L']*veh['m']*Ux[i]**2*K[i]) 
        aRHat = force2alpha(veh['FyRtable'], veh['alphaRearTable'], veh['a']/veh['L']*veh['m']*Ux[i]**2*K[i])    
        
        betaFFW.append( aRHat + veh['b']*K[i] )
        deltaFFW.append( K[i]*veh['L'] + aRHat - aFHat - veh['kLK']*veh['xLA']*betaFFW[i])
        
        [a, b, d, Cf, Cr] = getAffineModel(Ux[i], K[i], ts[i], veh, aFHat, aRHat, i);

        A[i] = a;
        B[i] = b;
        D[i] = d;
    
    print("DONE")
    
    return np.asarray(A), np.asarray(B), np.asarray(D), aF, aR, deltaFFW, betaFFW


def minimize_curvature(racing_line, Ux, inside, outside, safety_margin):
    
    print(" ---------- GET VEHICLE PARAMS -----------\n")
    
    ts = get_time(racing_line, Ux) # <<<<time

    N = len(racing_line); # number opt var	    
    
    s = get_s(racing_line)    
    
    s = s.reshape(1, len(s))

    veh = getVehicle('nonlinear', 'closest', N) # get vehicle params
    print("-Vehicle OK")
    
    print("\n ------------ GET ROAD WIDTH -------------- \n")
    width = get_width_with_rtree(inside, outside, racing_line) # get road width
    print("-Width OK")
    
    K   = get_K(racing_line) # get curvature
    
    print("\n --------------- GET SYSTEM DYNAMIC ------------- \n")
    
    [A, B, d, aF, aR, deltaFFW, betaFFW] = getAllSys(veh, Ux, K, ts) # System dynamic

    print("Dynamic OK")

    print("\n ------------- SETTING OPTIMIZATION PROBLEM ---------------- \n")

    opti = Opti(); # declaring optimization object
    
    # Symbolic variables

    X = opti.variable(5,N+1) #state trajectory X = [ e delta_psi r beta psi]'
         
    U = opti.variable(1,N) #control trajectory (steering)

    # Problem formulation
    
    J = sumsqr(diff(X[4,:])/(s)) + sumsqr(diff(U)); #cost function


    opti.minimize(J); # race in minimal time
    
    # Constraints
    perc = 0
    print("\nSetting Constraints...")
    for i in range(0,N):
        perc = print_perc(perc, N, i)
        opti.subject_to(X[:,i+1] - mtimes(A[i],X[:,i]) - mtimes(B[i],U[i]) - mtimes(d[i], ts[i]) == 0) # x_dot = A(t)*x + B(t)*delta + d(t) system dynamic
        opti.subject_to(opti.bounded(-width['inner_distance'][i] + safety_margin, X[0,i], width['outer_distance'][i] - safety_margin)) # road bounds
        

    opti.subject_to(X[0, 0] - X[0, -1] == 0) # continuity on e
    opti.subject_to(X[1, 0] + X[1, -1] == 0) # continuity on delta_psi
    opti.subject_to(X[2, 0] - X[2, -1] == 0) # continuity on yaw rate
    opti.subject_to(X[3, 0] - X[3, -1] == 0) # continuity on 
    opti.subject_to((X[0,-1] - X[0,-2])/s[0,-1] == (X[0,1] - X[0,0])/s[0,0])
    
    print("-Constraints OK")


    # Initial Guess

    # Solver
    popts = {"expand"  : True } # ipopt option
    sopts = {"max_iter" : 800} # ipopt option
    opti.solver('ipopt', popts, sopts) # setting solver
    sol = opti.solve() # solve proble

    opt_x = sol.value(X) # optimal states
    opt_y = sol.value(U) # optimal steering

    # Convert e and psi into global coordinates
    
    for i in range(0,N):
        

        racing_line[i,0] = racing_line[i,0] - opt_x[0,i]*np.cos(opt_x[4,i]);
        racing_line[i,1] = racing_line[i,1] - opt_x[0,i]*np.sin(opt_x[4,i]);
    
    
    return racing_line
    
    
    
    
    
    
