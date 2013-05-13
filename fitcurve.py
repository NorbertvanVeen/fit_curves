import scipy
import numpy as np
import math
import tables

from scipy.optimize import leastsq
from fit_form import coord_tcc, coin_stations
from scipy.stats import poisson, chi2
from scipy.special import gamma
from scipy.optimize import curve_fit
from numpy import nan, pi, sin, cos, tan, arctan2, arctan, arcsin, arccos, sqrt, histogram2d, vectorize, log10, average
from random import random

import random as random2
import scipy.optimize as sop
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

AGE_PARAMETER 			= .94
R_0 					= 40.
ALPHA					= 1.5
BETA					= 3.6
NUMBER_OF_PARTICLES     = [1e6,1e6]

def NKG_density1(distance, n_o_p):
	return n_o_p * normalising_factor(ALPHA, BETA) * (distance / R_0)**(AGE_PARAMETER-ALPHA) * (1+(distance/R_0))**(AGE_PARAMETER-BETA)

def makecoord(data, sel_coinc_ids, coords, t):
    c_index = data.root.coincidences.c_index
    waarde = c_index[t]
    readout = data.root.coincidences.observables.readCoordinates(waarde)
    readout_sort = sorted(readout, key=lambda o: o['station_id'])
    
    coinc_index = sel_coinc_ids.searchsorted(t)
    reconstructed_core = data.root.core_reconstructions.reconstructions.readCoordinates(coords, field='reconstructed_core_pos')
    x0, y0 = reconstructed_core[coinc_index] #showercore position
    print x0, y0
    
    lijst =[]
    gem_density =[]
    for event in readout_sort: # bepaling per event wat belangrijkste station is. i.e. the highest particle density
        lijst1 = (event['station_id'], event['n1'], event['n2'], event['n3'], event['n4'])
        gem = np.mean(lijst1[1:5])
        lijst.append(lijst1)
        gem_density.append(gem)
    
    print lijst
   
    #code on pc at home
    data2 = tables.openFile('C:/Python/station_coordinates.h5', 'r')
    distance2 =[]
    val_station = []
    for id, dens in enumerate(gem_density):
        x,y,z = data2.root.coordinates.stations[id][5]
        distance2.append(np.sqrt((x - x0)**2 + (y - y0)**2))
    data2.close()
    
    #code on laptop
    # dist_station =[]
    # cluster = data.root.core_reconstructions._v_attrs.cluster
    # for station in cluster.stations:
        # x, y, alpha = station.get_xyalpha_coordinates()
        # distance1 = np.sqrt((x - x1)**2 + (y - y1)**2)  
        # dist_station.append(distance1)
    
    
    return distance2, gem_density   #laptop: return dist_station, gem density

def fitcurve(x, y):
    popt, pcov = curve_fit(NKG_density1, x, y, [1300]) # fit met guess gebruikt werkt ook
        
    pylab.plot(NK_density1(range(0,600), *popt), '--')
    
    
def NKG_density(distance, n_o_p):
	r_rm = distance / R_0
	density = n_o_p * normalising_factor(ALPHA, BETA) * (r_rm)**(AGE_PARAMETER-ALPHA) * (1+r_rm)**(AGE_PARAMETER-BETA)
	return density

def normalising_factor(a, b):
	teller = gamma(b-AGE_PARAMETER)
	noemer = 2*pi * R_0**2 * gamma(AGE_PARAMETER-a+2) * gamma(a + b - (2*AGE_PARAMETER) - 2)
	return teller/noemer


def plot_NKG(Ne=[NUMBER_OF_PARTICLES[0]]):
    #plot function
    fig = plt.figure(figsize=(5.5,4.2))
    ax = fig.add_subplot(111, aspect='equal')
    #distance = np.arange(.1,SIM_r,.1)
    distance = np.arange(.1,200,.1)
    for n in Ne:
        density = [NKG_density(x,n_o_p=n) for x in distance]
        legend = 'log Ne = ' + str(round(math.log(n,10),1))
        ax.semilogy(distance, density, lw=2, label=legend)
        #ax.loglog(distance, density, label=legend)
    ax.grid(True)
    ax.legend()
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize='14')
    ax.set_xlabel('afstand tot showercore [m]',size='14')
    ax.set_ylabel('elektronendichtheid [m-2]',size='14')
    fig.show()

def plot_NKG2(distance2, dens):
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    distance = range(1,400)
    for macht in np.arange(4,8,1):
        n_o_p = 10**macht
        density = [NKG_density(x,n_o_p=n_o_p) for x in distance]
        legend = 'log Ne = ' + str(macht)
        ax.loglog(distance, density, label=legend)
    
    plt.scatter(distance2, dens, s = 5, color ='r', zorder=4)   
    ax.legend()
    ax.set_xlabel('afstand tot showercore [m]',size='14')
    ax.set_ylabel('elektronendichtheid [m-2]',size='14')
    ax.set_ylim(0.1,1e6)
    ax.set_xlim(0,600)
    ax.grid(True)
    fig.show()
    
    
if __name__ == '__main__':      
    n = 6
    tcc = 30
    
    with tables.openFile('C:/Python/results-mar-2013.h5', 'r') as data:
        coinc_ids, tccs, densities = coin_stations(data, n)
        sel_coinc_ids, result_coinc_ids, coords = coord_tcc(coinc_ids, tccs, data, tcc)
        print len(sel_coinc_ids)
        t = sel_coinc_ids[3] # zomaar een event met N==6 en Tcc>30
        dist2, dens = makecoord(data, sel_coinc_ids, coords, t)
        fitcurve(dist2, dens)
        plot_NKG2(dist2, dens)
        
        #plot_NKG2()
    