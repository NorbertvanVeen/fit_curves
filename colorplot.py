import numpy as np
import tables
import pylab

from math import *
from itertools import product, combinations
from scipy.optimize import leastsq
from progressbar import Bar, ETA, Percentage, ProgressBar

import matplotlib.pyplot as plt

def draw_plot_N6(data, n, tekst):
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
    
    
    print "Calculating plot ..."    
    widgets = [Percentage(), ' ', Bar(marker='#',left='[',right=']'), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets)
    
    print "number of events =", len(val_stat)
    
    for values in pbar(val_stat):
               
        col ='#fecc4f'
         

        plt.scatter(values[0][0], values[0][1], s =5, marker='.', color=col, zorder=1)
    
    #code to run on pc
    with tables.openFile('C:/Python/station_coordinates.h5', 'r') as data2:
        detector = data2.root.coordinates.stations
        for detect in detector:
            temp = [detect[1], detect[2], detect[3], detect[4]]
            station_id = [detect[0], detect[5]]
            for coord in temp:
                x,y = coord[0:2]
                plt.scatter(x , y, color ='r', zorder=4)
                x2,y2 = station_id[1][0]+14, station_id[1][1]-15
                ax.annotate(station_id[0], (x2,y2))
        data2.close()     
    
    
    #code for laptop
    # for station in cluster.stations:
        # for detector in station.detectors:
            # x, y, = detector.get_xy_coordinates()
            # plt.scatter(x , y, color ='r', zorder=4)
        # x, y, alpha = station.get_xyalpha_coordinates()
        # ax.annotate(station.station_id, (x, y))
        
    
    ax.set_title('reconstructed shower cores N == %s' %n)
    ax.set_xlabel('distance (m)')
    ax.set_ylabel('distance (m)')
    plt.xlim(-800, 800)
    plt.ylim(-800, 800)
    plt.show()
    return val_stat


def draw_plot_N5(data, n, tekst):
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
    k, a= 0,0
    
    print "Calculating plot ..."    
    widgets = [Percentage(), ' ', Bar(marker='#',left='[',right=']'), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets)
    
    print "number of events =", len(val_stat)
    
    for values in pbar(val_stat):
        
        if (values[1] ==True and values[2] ==True and values[3] ==True and values[4] ==True and values[6] ==True):
            col ='r'
            k +=1
               
        else:
            col ='#fecc4f'
            a +=1

        plt.scatter(values[0][0], values[0][1], s =5, marker='.', color=col, zorder=1)
    
    #code to run on pc
    with tables.openFile('C:/Python/station_coordinates.h5', 'r') as data2:
        detector = data2.root.coordinates.stations
        for detect in detector:
            temp = [detect[1], detect[2], detect[3], detect[4]]
            station_id = [detect[0], detect[5]]
            for coord in temp:
                x,y = coord[0:2]
                plt.scatter(x , y, color ='r', zorder=4)
                x2,y2 = station_id[1][0]+14, station_id[1][1]-15
                ax.annotate(station_id[0], (x2,y2))
        data2.close()     
    
    
    #code for laptop
    # for station in cluster.stations:
        # for detector in station.detectors:
            # x, y, = detector.get_xy_coordinates()
            # plt.scatter(x , y, color ='r', zorder=4)
        # x, y, alpha = station.get_xyalpha_coordinates()
        # ax.annotate(station.station_id, (x, y))
    
    print "kleur rood: station 501, 502, 503, 504, 506 = %s" %k
    
    print "goud: andere combinaties %s" %a 
    
    ax.set_title('reconstructed shower cores N == %s' %n)
    ax.set_xlabel('distance (m)')
    ax.set_ylabel('distance (m)')
    plt.xlim(-800, 800)
    plt.ylim(-800, 800)
    plt.show()
    return val_stat


def draw_plot_N(data, n, tekst):
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
    k, m, r, s, t, u, v, b, w, q, yi, z, h, a= 0,0,0,0,0,0,0,0,0,0,0,0,0,0 
    
    print "Calculating plot ..."    
    widgets = [Percentage(), ' ', Bar(marker='#',left='[',right=']'), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets)
    
    print "number of events =", len(val_stat)
    
    for values in pbar(val_stat):
        
        if (values[1] ==True and values[2] ==True and values[3] ==True):
            col ='r'
            k +=1
        elif (values[1] ==True and values[2] ==True and values[4] ==True):
            col ='c'
            m +=1
        elif (values[1] ==True and values[2] ==True and values[5] ==True):
            col ='m'
            r +=1
        elif (values[1] ==True and values[2] ==True and values[6] ==True):
            col ='g'
            s +=1
        elif (values[1] ==True and values[4] ==True and values[3] ==True):
            col ='#523f6d'
            t +=1
        elif (values[1] ==True and values[5] ==True and values[3] ==True):
            col ='#fdd300'
            u+=1
        elif (values[1] ==True and values[6] ==True and values[3] ==True):
            col ='k'
            v+=1
        elif (values[4] ==True and values[2] ==True and values[3] ==True):
            col ='#b4ffd6'
            b+=1
        elif (values[5] ==True and values[2] ==True and values[3] ==True):
            col ='#ae5013'
            w+=1
        elif (values[6] ==True and values[2] ==True and values[3] ==True):
            col = 'y'
            q+=1
        elif (values[3] ==True and values[4] ==True and values[5] ==True):
            col = '#ff1bba'
            yi+=1
        elif (values[4] ==True and values[5] ==True and values[6] ==True):
            col = '#840041'
            z+=1
        elif (values[3] ==True and values[4] ==True and values[6] ==True):
            col = '#f15a25'
            h+=1
        else:
            col ='#fecc4f'
            a +=1

        plt.scatter(values[0][0], values[0][1], s =5, marker='.', color=col, zorder=1)
    
    #code to run on pc
    with tables.openFile('C:/Python/station_coordinates.h5', 'r') as data2:
        detector = data2.root.coordinates.stations
        for detect in detector:
            temp = [detect[1], detect[2], detect[3], detect[4]]
            station_id = [detect[0], detect[5]]
            for coord in temp:
                x,y = coord[0:2]
                plt.scatter(x , y, color ='r', zorder=4)
                x2,y2 = station_id[1][0]+14, station_id[1][1]-15
                ax.annotate(station_id[0], (x2,y2))
        data2.close()     
    
    
    #code for laptop
    # for station in cluster.stations:
        # for detector in station.detectors:
            # x, y, = detector.get_xy_coordinates()
            # plt.scatter(x , y, color ='r', zorder=4)
        # x, y, alpha = station.get_xyalpha_coordinates()
        # ax.annotate(station.station_id, (x, y))
    
    print "kleur rood: station 501, 502, 503 = %s" %k
    print "kleur cyaan: station 501, 502, 504 = %s" %m
    print "kleur magenta: station 501, 502, 505 = %s" %r
    print "kleur groen: station 501, 502, 506 = %s" %s
    print "kleur paars: station 501, 503, 504 = %s" %t
    print "kleur donkergeel: station 501, 503, 505 = %s" %u
    print "kleur zwart: station 501, 503, 506 = %s" %v
    print "kleur mint: station 502, 503, 504 = %s" %b
    print "kleur donkerbruin: station 502, 503, 505 = %s" %w
    print "kleur geel: station 502, 503, 506 = %s" %q
    print "kleur fuchsia: station 503, 504, 505 = %s" %yi
    print "kleur bordeaux rood: station 504, 505, 506 = %s" %z
    print "kleur oranje: station 503, 504, 506 = %s" %h
    print "goud: andere combinaties %s" %a 
    
    ax.set_title('reconstructed shower cores N == %s' %n)
    ax.set_xlabel('distance (m)')
    ax.set_ylabel('distance (m)')
    ax.annotate('datafile =%s' %tekst, (-650, -650), color='k')
    plt.xlim(-800, 800)
    plt.ylim(-800, 800)
    plt.show()
    return val_stat
    
if __name__ == '__main__':  
    n  = 3
    tcc = 10
    tekst = 'results-alt-feb-2012.h5'
    with tables.openFile('C:/Python/%s' %tekst, 'r') as data:
        if (n <5):
            values = draw_plot_N(data,n, tekst)
        if (n==5):
            values = draw_plot_N5(data,n, tekst)
        elif(n==6):
            values = draw_plot_N6(data,n, tekst)