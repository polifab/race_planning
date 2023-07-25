
from math import pi
import numpy as np

def getVehicle(tireType, mapMatchType, N):
    
    nT = 250;
    N = nT    
    veh = {
        'a' : 1.7446,
        'b' : 1.1534,
        'm' : 1160.,
        'Cf' : 81277,
        'Cr' : 130470,
        'Iz' : 1260,
        'xLA' : 14.2,
        'kLK' : .0538,
        'muF' : 1.38,
        'muR' : 1.33,
        'g' : 9.81,
        'D' : .3638,
        'h' : .75,
        'alphaFlim' : 7*pi/180,
        'alphaRlim' : 5*pi/180
        }
    
    veh['L'] = veh['a'] + veh['b']
        
    veh['FzF'] = veh['m']*veh['b']*veh['g']/veh['L']
    veh['FzR'] = veh['m']*veh['a']*veh['g']/veh['L']
        
        
    veh['alphaFslide'] = abs( np.arctan2(3*veh['muF']*veh['m']*veh['b']/veh['L']*veh['g'],veh['Cf']) )
    veh['alphaRslide'] = abs( np.arctan2(3*veh['muR']*veh['m']*veh['a']/veh['L']*veh['g'],veh['Cr']) )
        
    #veh.alphaFrontTable:[-veh.alphaFslide:.001:veh.alphaFslide],   # vector of front alpha (rad)
    #veh.alphaRearTable :[-veh.alphaRslide:.001:veh.alphaRslide],   # vector of rear alpha (rad)
    veh['alphaFrontTable'] = np.linspace(-veh['alphaFslide'],veh['alphaFslide'], N) 
    #veh['alphaFrontTable'] = veh['alphaFrontTable'].reshape(1, len(veh['alphaFrontTable']))
    veh['alphaRearTable']  = np.linspace(-veh['alphaRslide'],veh['alphaRslide'], N)
    #veh['alphaRearTable']  = veh['alphaRearTable'].reshape(1, len(veh['alphaRearTable'])) 
    
    if tireType == 'linear':
        veh['FyFtable'] = -veh['Cf']*veh['alphaFrontTable'];
        veh['FyRtable'] = -veh['Cr']*veh['alphaRearTable'];
    elif tireType == 'nonlinear':
        veh['FyFtable'] = []
        veh['FyRtable'] = []
        for i in range(0,N):
            veh['FyFtable'].append(tireforces(veh['Cf'],veh['muF'],veh['muF'],veh['alphaFrontTable'][i],veh['FzF'], i));
            veh['FyRtable'].append(tireforces(veh['Cr'],veh['muR'], veh['muR'], veh['alphaRearTable'][i], veh['FzR'], i));
    else:
        error('Invalid Tire Type')
    
    veh['mapMatch'] = mapMatchType;
    veh['tireType'] = tireType;
    
    veh['brakeTimeDelay'] = .25; #seconds
    veh['rollResistance'] = 255; #Newtons
    veh['Kx'] = 3000; #speed tracking gain
    veh['powerLimit'] = 160000; #Watts
    
    return veh
    
    
def force2alpha(forceTable, alphaTable, Fdes): # mapping force to alpha

    if Fdes > max(forceTable):
        Fdes = max(forceTable) - 1;
    elif Fdes < min(forceTable):
        Fdes = min(forceTable) + 1;
    
    oldForce = forceTable[0]
 
    for i in range(0,len(forceTable)-1):
        
        if oldForce >= Fdes >= forceTable[i+1]:
            alpha = alphaTable[i]
            break
        else:
            oldForce = forceTable[i+1]
    
    return alpha



def tireforces(C,mu_peak,mu_slide,alpha,Fz, counter):
    # tireforces, calculate tire force for nonlinear tire model

    # Usage: calculate nonlinear tire curve, using Fiala brush model with no
    # longitudinal force.  Inputs are
    # C is the cornering stiffness per axle (N/rad)
    # mu_peak and mu_slide are no slip and slip friction coefficients
    # alpha = slip angle (rad)
    # Fz = normal load on that axle (N)

    alpha_slide = abs( np.arctan2(3*mu_peak*Fz,C) );
    
    
    if abs(alpha) < alpha_slide: # Not slide, use Fiala equations
        Fy = -C*np.tan(alpha)+(C)**(2./(3.*mu_peak*Fz))*(2-mu_slide/mu_peak)*abs(np.tan(alpha))*np.tan(alpha)-(C)**(3./(9.*(mu_peak)**2.*(Fz)**2.))*(np.tan(alpha))**3*(1-2*mu_slide/(3*mu_peak));
    else:    # Sliding on the surface
        Fy = -mu_slide*Fz*np.sign(alpha);
        
    return Fy



def getLocalStiffness(alpha, C, muP, muS, Fz, counter):

    delt = 1e-4;
    
    alpha2 = alpha + delt;
    alpha1 = alpha - delt;

    Fy1 = tireforces(C, muP, muS, alpha1, Fz, counter);
    Fy2 = tireforces(C, muS, muS, alpha2, Fz, counter);
    
    Cbar = (Fy2 - Fy1)/(alpha2 - alpha1);
    Cbar = abs(Cbar);
    
    return Cbar    
    
    
