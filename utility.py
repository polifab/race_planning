
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import math
from scipy import interpolate
from vehicle import *
from boostrtrees import RTree
from scipy.linalg import expm
from scipy.signal import cont2discrete

def get_time(racing_line, speed_profile):
    return np.sqrt(np.gradient(racing_line[:,0])**2 + np.gradient(racing_line[:,1])**2)/speed_profile

def interp(x, rate):
    
    N = len(x)
    
    tck = interpolate.splrep(np.arange(N), x, s=0)
    
    xx = np.linspace(0, N, N*rate) 

    return interpolate.splev(xx, tck, der=0)

def new_interp_line(line, rate):
    
    new_x = interp(line[:,0], int(rate))
    new_y = interp(line[:,1], int(rate))
    
    new_x = new_x.reshape(len(new_x), 1)
    new_y = new_y.reshape(len(new_y), 1)
    
    new_line = np.concatenate((new_x, new_y), axis=1) 
    
    return new_line 



def downsample(x,rate):
    
    new_x = np.zeros(int(len(x)*rate))
    
    for i in range(0,len(new_x)):
        new_x[i] = x[int(i/rate)]

    return new_x

    
def new_down_line(line, rate):
    
    new_x = downsample(line[:,0], rate)
    new_y = downsample(line[:,1], rate)

    new_x = new_x.reshape(len(new_x), 1)
    new_y = new_y.reshape(len(new_y), 1)
    
    new_line = np.concatenate((new_x, new_y), axis=1) 
    
    return new_line


def preprocessing(racing_line, rate = 6, k_reg = 3, s_reg = 10, stepsize_prep = 1.0, stepsize_reg = 3.0):
        
        track_cl = np.vstack((racing_line, racing_line[0]))
        no_points_track_cl = track_cl.shape[0]

        # calculate element lengths (euclidian distance)
        el_lengths_cl = np.sqrt(np.sum(np.power(np.diff(track_cl[:, :2], axis=0), 2), axis=1))

        # sum up total distance (from start) to every element
        dists_cum_cl = np.cumsum(el_lengths_cl)
        dists_cum_cl = np.insert(dists_cum_cl, 0, 0.0)

        # calculate desired lenghts depending on specified stepsize (+1 because last element is included)
        no_points_interp_cl = math.ceil(dists_cum_cl[-1] / stepsize_prep) + 1
        dists_interp_cl = np.linspace(0.0, dists_cum_cl[-1], no_points_interp_cl)

        # interpolate closed track points
        track_interp_cl = np.zeros((no_points_interp_cl, 2))
        track_interp_cl[:, 0] = np.interp(dists_interp_cl, dists_cum_cl, track_cl[:, 0])
        track_interp_cl[:, 1] = np.interp(dists_interp_cl, dists_cum_cl, track_cl[:, 1])
        #track_interp_cl[:, 2] = np.interp(dists_interp_cl, dists_cum_cl, track_cl[:, 2])
        #track_interp_cl[:, 3] = np.interp(dists_interp_cl, dists_cum_cl, track_cl[:, 3])

        # ------------------------------------------------------------------------------------------------------------------
        # SPLINE APPROXIMATION / PATH SMOOTHING ----------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------------------

        # find B spline representation of the inserted path and smooth it in this process
        # (tck_cl: tuple (vector of knots, the B-spline coefficients, and the degree of the spline))
        tck_cl, t_glob_cl = interpolate.splprep([track_interp_cl[:, 0], track_interp_cl[:, 1]],
                                                k=k_reg,
                                                s=s_reg,
                                                per=1)[:2]

        # calculate total length of smooth approximating spline based on euclidian distance with points at every 0.25m
        no_points_lencalc_cl = math.ceil(dists_cum_cl[-1]) * 4
        path_smoothed_tmp = np.array(interpolate.splev(np.linspace(0.0, 1.0, no_points_lencalc_cl), tck_cl)).T
        len_path_smoothed_tmp = np.sum(np.sqrt(np.sum(np.power(np.diff(path_smoothed_tmp, axis=0), 2), axis=1)))

        # get smoothed path
        no_points_reg_cl = math.ceil(len_path_smoothed_tmp / stepsize_reg) + 1
        racing_line = np.array(interpolate.splev(np.linspace(0.0, 1.0, no_points_reg_cl), tck_cl)).T[:-1]
        
        racing_line = new_interp_line(racing_line, rate)
        
        return racing_line
        
        
        
        
        
        
        
        
def print_perc(perc, N, i):
    
    if round(i*100./N) != perc  and int(round(i*100./N)) != 100:
        perc = round(i*100./N)
        b = str(">> "+str(int(perc)))+"%"  
        print (b, end="\r")
    elif round(i*100./N) != perc and int(round(i*100./N)) == 100:
        print(">> 100%\n")
        perc = round(i*100./N)

    return perc


def dist(a,b):
    return np.sqrt((b[1] - a[1])**2 + (b[0] - a[0])**2)

def find_nearest(out_i, inn, counter):
    
    minimum = 10e9
    count = 0
        
    for i in range(counter - 300, counter + 300):
#        if i < 0:
#            i = 0
        if abs(i) >= len(inn):
            i = np.sign(i)*len(inn) -1
        
        if dist(out_i, inn[i]) <= minimum:
            minimum = dist(out_i, inn[i])
            min_point = inn[i]
            
    return min_point, minimum

def get_width(inner_bound, outer_bound, racing_line):

    vec_len = len(racing_line)     
    out_width = np.zeros(vec_len)
    in_width = np.zeros(vec_len)
    perc = 0
    print("-Computing road width...")
    for i in range(0, vec_len):
        perc = print_perc(perc, vec_len, i)           
        outbound_corr = len(outer_bound)/float(vec_len)
        [out_min_point, out_width[i]] = find_nearest(racing_line[i], outer_bound, int(np.floor(outbound_corr * i)))

        inbound_corr = len(inner_bound)/float(vec_len)
        [in_min_point, in_width[i]] = find_nearest(racing_line[i], inner_bound, int(np.floor(inbound_corr * i)))			

    #print("DONE")
		    
    width = {'outer_distance': out_width, 'inner_distance': in_width} 

    return width	


def get_width_with_rtree(inner_bound, outer_bound, racing_line):
    vec_len = len(racing_line)     
    out_width = np.zeros(vec_len)
    in_width = np.zeros(vec_len)
    perc = 0
    print("-Computing road width...")
    print(" -> step_1: creating rtree structures...\n")
    # inner_bound
    inner_rtree = RTree()
    length = len(inner_bound) 
    for i in range(length):
        # insert point (x-value,y-value, index inside array)
        inner_rtree.insert_point(inner_bound[i][0],inner_bound[i][1], i)
    print("         inner_bound.....OK\n")
    # outer_bound
    outer_rtree = RTree()
    length = len(outer_bound) 
    for i in range(length):
        # insert point (x-value,y-value, index inside array)
        outer_rtree.insert_point(outer_bound[i][0],outer_bound[i][1], i)
    print("         outer_bound.....OK\n")
    
    print(" -> step_2: find min distance...")
    for i in range(0, vec_len):
        perc = print_perc(perc, vec_len, i) 
        out_width[i]=outer_rtree.min_distance(racing_line[i][0],racing_line[i][1])
        in_width[i]=inner_rtree.min_distance(racing_line[i][0],racing_line[i][1])

    width = {'outer_distance': out_width, 'inner_distance': in_width}
    return width
    
    
def get_s(racing_line):
    return (np.sqrt(np.gradient(racing_line[:,0])**2 + np.gradient(racing_line[:,1])**2))

def get_K(racing_line):
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


def get_yaw(racing_line):
    diffX = np.gradient(racing_line[:,0]); 
    diffY = np.gradient(racing_line[:,1]);
    psi_temp = np.arctan2(diffY, diffX);
    return psi_temp
