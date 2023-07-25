# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 15:08:13 2020

@author: Fabio
"""
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def gLongMax(gLat_current,gLat_Max, gLong_Max):

	a = 1 / (gLat_Max**2)
	b = 1 / (gLong_Max**2)

	max_gLong = abs(np.sqrt(abs((1 - a*(gLat_current**2))/b)))

	return max_gLong	


def get_curvature(racing_line):
	x = racing_line[:,0]
	y = racing_line[:,1]
	dx  = np.gradient(x)
	ddx = np.gradient(dx)
	dy  = np.gradient(y)
	ddy = np.gradient(dy)
	num   = dx * ddy - ddx * dy
	denom = dx * dx + dy * dy	
	denom = denom**(3.0/2.0)
	curv = num / denom

	return curv


def get_SpeedProfile(Racing_Line, glat_max, glong_acc_max, glong_dec_max, vMax):

        Racing_Line_K = get_curvature(Racing_Line);   
    
        mu_x = glat_max;
        gLong_acc = glong_acc_max*9.8;
        gLong_dec = glong_dec_max*9.8;
        max_speed = vMax/3.6;
        max_acc = mu_x*9.8;
        #g_correct = 9.645; # valore di g corretto al ribasso per limite massimo sulle gLong
        #gLong_acc = options.gLong_acc #*g_correct;
        #gLong_dec = options.gLong_dec #*g_correct;

        #lat_coeff_acc = 0.1;
        #lat_coeff_dec = 0.1;
        X = Racing_Line[:,0];
        Y = Racing_Line[:,1];

        diff_X = np.gradient(X);
        diff_Y = np.gradient(Y);
        s = np.sqrt((diff_X)**2 + (diff_Y)**2);
        #len_s = len(s);
        #S = np.cumsum(s);        
        
        
        #Psi = np.arctan2(diff_Y, diff_X);
        #Psi_unw = np.unwrap(Psi);  
        K = Racing_Line_K;
           
        #calcolo vmax
        v1square = (abs(K)**(-1))*max_acc;
        v1 = np.sqrt(v1square);
                     
        for i in range(0,len(v1)):
            if v1[i] > max_speed:
                v1[i] = max_speed
           
        # forward sweep (acceleration limit)
        v2 = v1;
        actual_velocity = v2[0];
        actual_aN = (actual_velocity**2)*abs(K[0]);
        actual_aT = gLongMax(actual_aN, max_acc, gLong_acc);
        #curvature = Racing_Line_K;
        #percent_K = max(abs(curvature))*0.10;
        #max_acc_old = max_acc;

        for j in range(1,max(v2.shape)):
            
            # per avere accelerazione costante in curva
            #if(abs(curvature(j)) >= percent_K)
            #     max_acc = 0;
            #  else
            #      max_acc = max_acc_old;
            #  end
            #--------------------------
            
            next_velocity = v1[j];
            
            #v2(j) = min(next_velocity, actual_velocity.*realsqrt(1 + 2.*s(j)*actual_aT./(actual_velocity.^2)));

            v2[j] = min(next_velocity, actual_velocity*np.sqrt(1 + 2*s[j]*actual_aT/float(actual_velocity**2)));
         
            actual_velocity = v2[j];
            actual_aN = (actual_velocity**2)*abs(K[j]);
            actual_aT = gLongMax(actual_aN, max_acc, gLong_acc);
        

        # BACKWARD SWEEP (DECELERATION LIMIT)
        v3 = v2;
        next_velocity = v3[len(v3)-1];
        
        al = np.array([])
        at = np.array([])
        
        next_aN = (next_velocity**2)*abs(K[len(K)-1]);
        next_aT = gLongMax(next_aN, max_acc, gLong_dec); #prendo il limite sull'ellisse

        al = np.append(al,next_aN)    
        at = np.append(at,next_aT)
        for j in reversed(range(1, max(v3.shape))):
        
            #  per avere accelerazione costante in curva
            #  if(abs(curvature(j)) >= percent_K)
            #      max_acc = 0.1;
            #  else
            #      max_acc = max_acc_old;
            #  end
            #
            #--------------------------
            actual_velocity = v2[j-1];
            v3[j-1] = min(actual_velocity, next_velocity*np.sqrt(1 + 2*s[j-1]*next_aT/float(next_velocity**2)));
            next_velocity = v3[j-1];
            next_aN = (next_velocity**2)*abs(K[j]);
            next_aT = gLongMax(next_aN, max_acc, gLong_dec);
        
            al = np.append(al,next_aN)
            at = np.append(at,next_aT)
        
        return v3

def get_SpeedProfile_rect(Racing_Line, glat_max, glong_acc_max, glong_dec_max, vMax):

        Racing_Line_K = get_curvature(Racing_Line);   
    
        mu_x = glat_max;
        gLong_acc = glong_acc_max*9.8;
        gLong_dec = glong_dec_max*9.8;
        max_speed = vMax/3.6;
        max_acc = mu_x*9.8;
        #g_correct = 9.645; # valore di g corretto al ribasso per limite massimo sulle gLong
        #gLong_acc = options.gLong_acc #*g_correct;
        #gLong_dec = options.gLong_dec #*g_correct;

        #lat_coeff_acc = 0.1;
        #lat_coeff_dec = 0.1;
        X = Racing_Line[:,0];
        Y = Racing_Line[:,1];

        diff_X = np.gradient(X);
        diff_Y = np.gradient(Y);
        s = np.sqrt((diff_X)**2 + (diff_Y)**2);
        #len_s = len(s);
        #S = np.cumsum(s);        
        
        
        #Psi = np.arctan2(diff_Y, diff_X);
        #Psi_unw = np.unwrap(Psi);  
        K = Racing_Line_K;
           
        #calcolo vmax
        v1square = (abs(K)**(-1))*max_acc;
        v1 = np.sqrt(v1square);
                     
        for i in range(0,len(v1)):
            if v1[i] > max_speed:
                v1[i] = max_speed
           
        # forward sweep (acceleration limit)
        v2 = v1;
        actual_velocity = v2[0];
        actual_aN = (actual_velocity**2)*abs(K[0]);
        actual_aT = gLong_acc
        #curvature = Racing_Line_K;
        #percent_K = max(abs(curvature))*0.10;
        #max_acc_old = max_acc;

        for j in range(1,max(v2.shape)):
            
            # per avere accelerazione costante in curva
            #if(abs(curvature(j)) >= percent_K)
            #     max_acc = 0;
            #  else
            #      max_acc = max_acc_old;
            #  end
            #--------------------------
            
            next_velocity = v1[j];
            
            #v2(j) = min(next_velocity, actual_velocity.*realsqrt(1 + 2.*s(j)*actual_aT./(actual_velocity.^2)));

            v2[j] = min(next_velocity, actual_velocity*np.sqrt(1 + 2*s[j]*actual_aT/float(actual_velocity**2)));
         
            actual_velocity = v2[j];
            actual_aN = (actual_velocity**2)*abs(K[j]);
            actual_aT = gLong_acc
        

        # BACKWARD SWEEP (DECELERATION LIMIT)
        v3 = v2;
        next_velocity = v3[len(v3)-1];
        
        al = np.array([])
        at = np.array([])
        
        next_aN = (next_velocity**2)*abs(K[len(K)-1]);
        next_aT = gLong_dec #prendo il limite

        al = np.append(al,next_aN)    
        at = np.append(at,next_aT)
        for j in reversed(range(1, max(v3.shape))):
        
            #  per avere accelerazione costante in curva
            #  if(abs(curvature(j)) >= percent_K)
            #      max_acc = 0.1;
            #  else
            #      max_acc = max_acc_old;
            #  end
            #
            #--------------------------
            actual_velocity = v2[j-1];
            v3[j-1] = min(actual_velocity, next_velocity*np.sqrt(1 + 2*s[j-1]*next_aT/float(next_velocity**2)));
            next_velocity = v3[j-1];
            next_aN = (next_velocity**2)*abs(K[j]);
            next_aT = gLong_dec
        
            al = np.append(al,next_aN)

            at = np.append(at,next_aT)
        
        return v3

if __name__ == "__main__":
    
    print("Script under construnction...")
