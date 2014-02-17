import numpy as np
import tables
import pylab

from math import *
from itertools import product, combinations
from scipy.optimize import leastsq

import matplotlib.pyplot as plt





# aangepast op 21 maart, nu bruikbaar voor events die een coincididentie van minder dan 6 stations hebben

def coin_stations(data, n):
   
    coinc_ids = data.root.coincidences.coincidences.readWhere('N == %d' %n, field='id')
    #id = coincidences['id'] #where coincidences of 3 stations are placed
    c_index = data.root.coincidences.c_index
    observ = data.root.coincidences.observables
    
    n = []
    tccs= []
    plate_dens =[]
    for coinc_id in coinc_ids:
        event_ids = c_index[coinc_id]
        events = observ.readCoordinates(event_ids)
        densities = []
        for event in events:
            event_densities =[event[u] for u in 'n1', 'n2', 'n3', 'n4']
            densities.extend(event_densities)        
        plate_dens.extend(densities)
        n = len(densities)            #bepaling van het aantal platen wat mee doet aan de bepaling
        n_mean = np.mean(densities)   #bepaling TCC volgens Jos Steijger, gemiddelde 
        S_var = np.var(densities)     #bepaling TCC volgens Jos Steijger, variantie        
        Ts = S_var / n_mean
        Tcc = (n - 1) * Ts
        tccs.append(Tcc)       
        
    return coinc_ids, np.array(tccs), plate_dens   

    
def coord_tcc(coinc_ids, tccs, data, tcc): 
    # coinc_ids = data.root.coincidences.coincidences.readWhere('N >= 3', field= 'id')  # plek van de coincidenties vinden waar de TCC waarde hoger is dan 30
    sel_coinc_ids = coinc_ids.compress(tccs > tcc)  # voor bepaling voor TCC waarden hoger of lager dan vooraf ingestelde waarde
    # print len(sel_coinc_ids)
    
    N = len(data.root.core_reconstructions.reconstructions)
    rec_ids = np.arange(N)
    rec_coinc_ids = data.root.core_reconstructions.reconstructions.col('coinc_id')   # truc om de plekken van de coincidenties te vinden.

    coords = rec_ids.compress([True if id in sel_coinc_ids else False for id in rec_coinc_ids]) #coordinaten van de plekken

    result_coinc_ids = data.root.core_reconstructions.reconstructions.readCoordinates(coords, field='coinc_id')
    print "lengte result_coinc_ids:", len(result_coinc_ids)
    print "lengte sel_coinc_ids:", len(sel_coinc_ids)
    #result coinc ids should be equal to sel_coinc_ids
    return sel_coinc_ids, result_coinc_ids, coords

    
def draw_stations_tcc(data, coords, n, tcc):
    cluster = data.root.core_reconstructions._v_attrs.cluster
    fig = plt.figure()
    ax = fig.add_subplot(111)
     
    core_tcc = data.root.core_reconstructions.reconstructions.readCoordinates(coords, field='reconstructed_core_pos') #levert ....showers
       
    list_core =[]
    for i in core_tcc:
        list_c =[]
        if (abs(i[0])<1000) & (abs(i[1]) <1000):
            list_c = [i[0], i[1]]
            list_core.append(list_c)
    list_core = np.array(list_core)  #... cores die goed zijn.
    plt.plot(list_core[:,0], list_core[:,1], ',', zorder=1) #geeft een plot van de core posities
    
    # for station in cluster.stations:
        # x,y, alpha = station.get_xyalpha_coordinates()
        
    for station in cluster.stations:
        for detector in station.detectors:
            x, y, = detector.get_xy_coordinates()
            plt.scatter(x , y, color ='r', zorder=4)
        x, y, alpha = station.get_xyalpha_coordinates()
        ax.annotate(station.station_id, (x, y))
    
    ax.set_title('reconstructed shower cores N == %s' %n)
    ax.set_xlabel('distance (m)')
    ax.set_ylabel('distance (m)')
    ax.annotate('number of coincidences %s' %len(coords), (-700, -700))
    ax.annotate('tcc-value >= %s' %tcc, (-650, -650))
    plt.xlim(-800, 800)
    plt.ylim(-800, 800)
    plt.show()
    


def plot_core(data, sel_coinc_ids):
    
    cluster = data.root.core_reconstructions._v_attrs.cluster
    c_index = data.root.coincidences.c_index
    observ = data.root.coincidences.observables
    #coinc_ids = data.root.coincidences.coincidences.readWhere('N==3', field='id')
    #sel_coinc_ids = coinc_ids.compress(tccs > 30)
    #waarde = c_index[sel_coinc_ids]  # find the place in observables.
    #print sel_coinc_ids
    sel_coinc_ids1 = sel_coinc_ids[1:4]
    for i in sel_coinc_ids1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        coinc_index = sel_coinc_ids.searchsorted(i)
        print coinc_index
        print i
        print sel_coinc_ids.shape

        print "event nummer =", i
        waarde = c_index[i]
        readout = data.root.coincidences.observables.readCoordinates(waarde)
        readout_sort = sorted(readout, key=lambda o: o['station_id']) #sorteert op station id
        #print readout_sort
        lijst =[]
        gem_density =[]
        station_nr =[]
        for event in readout_sort: # bepaling per event wat belangrijkste station is. i.e. the highest particle density
            lijst1 = (event['station_id'], event['n1'], event['n2'], event['n3'], event['n4'])
            gem = np.mean(lijst1[1:5])
            nr = event['station_id']
            lijst.append(lijst1)
            gem_density.append(gem)
            station_nr.append(nr)
        print "lijst met densities + station",lijst
        
        station_nr = np.array(station_nr)
        print "station nummers die meedoen aan het event:", station_nr
        tot_dens=[]
        for rij in lijst:
            tot = np.sum(rij[1:5])
            tot_dens.append(tot) #geeft een lijst met een optelling van de particle densities
        
        index_station = np.array(tot_dens.index(max(tot_dens))) # Max waarde uit de lijst is ook waarde van station id met deze max waarde
        val_maxstation = lijst[index_station]
        val_maxstation1 = val_maxstation[1:5]
        
        print "lijst van totale densities", tot_dens
        print "index_station", index_station
        print "val max", val_maxstation1 
       
       # code om de grootste waarde te krijgen uit de lijst. Bepaal r ten opzichte van die coordinaten.
        # value_id = [x for y, x in zip(tot_dens, lijst) if y == max(tot_dens)]
        # print "value id is", value_id
        # id_stat = value_id[0][0]
        # bovenstaande code is niet nodig omdat we de positie van de core nu pakken uit de lijst
        #ax.annotate(station.station_id, (x, y))
        #x0, y0, alpha0 = cluster.stations[id_stat].get_xyalpha_coordinates() #place the core here closest to the station with highest plate density.
        
        reconstructed_core = data.root.core_reconstructions.reconstructions.readCoordinates(coords, field='reconstructed_core_pos')
        # rec = data.root.core_reconstructions.reconstructions.readCoordinates(coords, field='coinc_id')
        # print reconstructed_core.shape
        x0, y0 = reconstructed_core[coinc_index]        # coordinaten van reconstructed shower core element is het element van de gekozen sel_coinc_ids1
        # x0 = 179
        # y0 = 118
        print 80*'-'
        # print i
        # print coinc_index
        # print rec[coinc_index]
       # print 80*'-'
        
        dist_detect = []
        station = cluster.stations[index_station]
        for detector in station.detectors:
            x, y = detector.get_xy_coordinates()
            distance2 = np.sqrt((x - x0)**2 + (y - y0)**2)
            dist_detect.append(distance2)
                        
        print "distance detectors", dist_detect
                
               
        dist_station =[]
        for station in [cluster.stations[i] for i in station_nr]:
            x, y, alpha = station.get_xyalpha_coordinates()
            print "stationnr en coordinaten =", i , x, y
            distance1 = np.sqrt((x - x0)**2 + (y - y0)**2)  
            dist_station.append(distance1)
        tot =[x0, y0]       
        gem_density = sorted(gem_density, reverse=True)
        dist_station = sorted(dist_station)
        
        print 80*'-'
        print "lijst met distance station", dist_station
        print "lijst met gem density",gem_density
        print "coordinaten core x, y", x0, y0
         
        #for i in range(0,6):
        #    plt.scatter(dist_station[i], gem_density[i], zorder=4)
        
        #plt.text(200, 7,'uitgerekend met Visuele showercore op %s' %tot)
       
        plt.scatter(dist_detect, val_maxstation1, color ='r', zorder=2)
        plt.scatter(dist_station, gem_density)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_title('distance versus particle density m-2 event %s' % i)
        ax.set_xlabel('R (m)')
        ax.set_ylabel('particle density ')
        
        plt.xlim(0, 1400)
        # plt.ylim(0, 11)
        plt.show()      
        
            
            
if __name__ == '__main__':  
    n  = 4
    tcc = 10
    data = tables.openFile('C:/Python/results-mar-2013.h5', 'r')
    coinc_ids, tccs, densities = coin_stations(data, n)
    sel_coinc_ids, result_coinc_ids, coords = coord_tcc(coinc_ids, tccs, data, tcc)
    draw_stations_tcc(data, coords, n, tcc)
    print "coincidenties %s voudig" %n, len(sel_coinc_ids)
    data.close()
    #plot_core(data, sel_coinc_ids)
