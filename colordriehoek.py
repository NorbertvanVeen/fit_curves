import numpy as np
import tables
import pylab

from math import *
from itertools import product, combinations
from scipy.optimize import leastsq
from progressbar import Bar, ETA, Percentage, ProgressBar

import matplotlib.pyplot as plt


def draw_plot_N(data, n, stat1, stat2, stat3, tekst):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cluster = data.root.core_reconstructions._v_attrs.cluster
    list_coinc =data.root.core_reconstructions.reconstructions
    list = list_coinc.getWhereList('N==n')
    readout = list_coinc.readCoordinates(list)
    
    core_pos =[]
    val_stat =[]
    for l, row in enumerate(readout):
        temp1 = []
        temp1 = row[3], row[5], row[6], row[7], row[8], row[9], row[10]
        val_stat.append(temp1)
    
    k,a  = 0, 0
    
    print "Calculating plot ..."    
    widgets = [Percentage(), ' ', Bar(marker='#',left='[',right=']'), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets)
        
    for values in pbar(val_stat):
        if (values[stat1] ==True and values[stat2] ==True and values[stat3] ==True):
            col ='b'
            k +=1
            
            plt.scatter(values[0][0], values[0][1], s =5, marker='.', color=col, zorder=1)
        else:
            a+=1
    #code to run on pc
    with tables.openFile('C:/Python/station_coordinates.h5', 'r') as data2:
        detector = data2.root.coordinates.stations
        for detect in detector:
            temp = [detect[1], detect[2], detect[3], detect[4]]
            station_id = [detect[0], detect[5]]
            for coord in temp:
                x,y = coord[0:2]
                plt.scatter(x , y, color ='r', zorder=4)
                x2,y2 = station_id[1][0], station_id[1][1]
                ax.annotate(station_id[0], (x2,y2))
        data2.close()        
                
    # for station in cluster.stations:
        # for detector in station.detectors:
            # x, y, = detector.get_xy_coordinates()
            # plt.scatter(x , y, color ='r', zorder=4)
        # x, y, alpha = station.get_xyalpha_coordinates()
        # ax.annotate(station.station_id, (x, y))
    
    station = [stat1, stat2, stat3]
    print "gekozen stations", stat1, stat2, stat3
    print "number of showers: = %s" %k
    
    print "andere combinaties %s" %a
    
    ax.annotate('stations =%s' %station, (-700, -650))
    ax.annotate('aantal reconstructies %s' %k, (-700, -700))
    ax.annotate('datafile =%s' %tekst, (-700, -600), color='k')
    ax.set_title('reconstructed shower cores N == %s' %n)
    
    ax.set_xlabel('distance (m)')
    ax.set_ylabel('distance (m)')
    plt.xlim(-800, 800)
    plt.ylim(-800, 800)
    plt.show()
    return val_stat
    
if __name__ == '__main__':  
    n  = 3
    stat1 = 1
    stat2 = 2
    stat3 = 6
    tekst = 'results-alt-feb-2012.h5'
    tcc = 10
    
    with tables.openFile('C:/Python/%s' %tekst, 'r') as data:
        values = draw_plot_N(data, n, stat1, stat2, stat3, tekst)
        
    