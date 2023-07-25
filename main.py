
import json
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import math
from scipy import interpolate
import datetime

from three_steps_SpeedProfile import get_SpeedProfile
from three_steps_SpeedProfile import get_SpeedProfile_rect
from fast_racing_line import minimize_curvature

from utility import *

if __name__ == '__main__':
    
    
    if  len(sys.argv) != 9   :
        print ('***INVALID n# of args !***')
        print ('Specify: track_name')
        print ('Closing')
        exit(0)

    else:    
        
        track_name  = str(sys.argv[1])
        
        track_name  = str(sys.argv[1])
        name        = str(sys.argv[2])
        vMax        = float(sys.argv[3])      # km/h
        gAcc        = float(sys.argv[4])      # g
        gDec        = float(sys.argv[5])      # g
        gLat        = float(sys.argv[6])      # g
        safety_margin   = float(sys.argv[7])  # m
        
        output_destination = str(sys.argv[8])
        
        
        start_time = time.time()
        
        print("\n\n******************************************************************************************************************************************************")
        
        print(" ______   ______     ______     ______      ______     ______     __   __     ______     ______     ______     ______   __     ______     __   __")    
        print("/\  ___\ /\  __ \   /\  ___\   /\__  _\    /\  ___\   /\  ___\   /\ \-.\ \   /\  ___\   /\  == \   /\  __ \   /\__  _\ /\ \   /\  __ \   /\ \-.\ \ ")  
        print("\ \  __\ \ \  __ \  \ \___  \  \/_/\ \/    \ \ \__ \  \ \  __\   \ \ \-.\ \  \ \  __\   \ \  __<   \ \  __ \  \/_/\ \/ \ \ \  \ \ \/\ \  \ \ \-.\ \ ") 
        print(" \ \_\    \ \_\ \_\  \/\_____\    \ \_\     \ \_____\  \ \_____\  \ \_\  \_\  \ \_____\  \ \_\ \_\  \ \_\ \_\    \ \_\  \ \_\  \ \_____\  \ \_\  \_\ ")
        print("  \/_/     \/_/\/_/   \/_____/     \/_/      \/_____/   \/_____/   \/_/ \/_/   \/_____/   \/_/ /_/   \/_/\/_/     \/_/   \/_/   \/_____/   \/_/ \/_/ ")
        
        print("\n******************************************************************************************************************************************************\n\n")


        print("----- Read track from json file --------\n")
        with open(track_name, 'r') as json_rac:
            rac_dict = json.load(json_rac)
        json_rac.close()
        
        racing_line = np.array(rac_dict['Centre'])
        old_racing = np.copy(racing_line)
        inside = np.array(rac_dict['Inside'])
        outside = np.array(rac_dict['Outside'])
        
        print("-JSON OK\n")

        inside = new_down_line(inside, 2)
        outside = new_down_line(outside, 2)
        
        inside = preprocessing(inside)
        outside = preprocessing(outside)
        
        #---------------- PREPROCESSING LINES ------------------------------- 
        
        racing_line = preprocessing(racing_line)
        
     
        
        print("-PREPROCESSING OK\n")
        # ..........................................................
        
        
        max_iteration   = 10
        n_iterations    = 0
        
        
        minimum = 10e9
        best_lap = 0
        delta_t = 10e9
        
        while n_iterations < max_iteration:
            
            
            print("************** STARTING WITH " + str(n_iterations+1) + "° LOOP OF OPTIMIZATION ****************\n")

            # new speed profile

            speed_profile   =   get_SpeedProfile_rect(racing_line, gLat, gAcc, gDec, vMax)

            t   =   sum(get_time(racing_line, speed_profile))
            
            # -----------------

            # new racing_line
            
            racing_line     = minimize_curvature(racing_line, speed_profile, inside, outside, safety_margin)
            speed_profile   = get_SpeedProfile_rect(racing_line, gLat, gAcc, gDec, vMax)
            
            # ---------------
            
            t_old   = t
            t       = sum(get_time(racing_line, speed_profile))
            delta_t = t_old - t
            K = get_K(racing_line)
            time_elapsed = time.time() - start_time
            
            print("- TIME = " + str(t))
            print("- K = " + str(sum(np.abs(K))))
            print("- DELTA T = " + str(delta_t) + " (positive means an improvement)")       
            print("- TIME ELAPSED == " + str(time_elapsed))

            if t < minimum:
                best_racing_line = racing_line
                best_speed_profile = speed_profile
                minimum = t
                best_lap = n_iterations + 1
                
            n_iterations += 1

        print("------------------------- OPTIMIZATION FINISHED ---------------------------------------")
        print("BEST LAP IS: " + str(best_lap) + " WITH A LAP TIME OF " + str(minimum) + " s")
        
        fig = plt.figure()
        plt.plot(outside[:,0], outside[:,1])        
        plt.plot(inside[:,0], inside[:,1])
        plt.plot(best_racing_line[:,0], best_racing_line[:,1])
        plt.title("Racing Line")
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.show()
        fig.savefig('racing_line_'+str(n_iterations+1)+'.png')

        plt.plot(speed_profile)
        plt.title("Speed Profile")
        plt.xlabel('point (k)')
        plt.ylabel('speed (m/s)')
        plt.show()

        final_time = time.time() - start_time
        print("TOTAL TIME ELAPSED IS + " + str(final_time) + "\n\n")

        data = {
            "Date" : str(datetime.datetime.now()),
            "Params" : [vMax, gAcc, gDec, gLat, safety_margin],
            "r_l" : np.ndarray.tolist(best_racing_line),
            "s_p" : np.ndarray.tolist(best_speed_profile)
        }

        with open(output_destination + "//full_fst_" + name + "_m" + str(safety_margin) + ".json", 'w') as outfile:
            json.dump(data, outfile)
        
        
        
