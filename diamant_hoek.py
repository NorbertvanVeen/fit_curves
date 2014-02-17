import os

import tables
import pylab
import scipy.ndimage
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from numpy import arctan2, cos, sin, arcsin, isnan, pi, linspace
from scipy.optimize import curve_fit
from pylab import histogram

from hisparc.analysis.traces import get_traces
from sapphire.analysis.process_events import ProcessEvents, ProcessIndexedEvents
from sapphire.analysis.direction_reconstruction import DirectionReconstruction

import numpy as np

STATION_ID = [501, 502, 503, 504, 505, 506, 508]
DATAPATH = "/Users/cgvanveen/HiSPARC/data"

# dit programma doet:
# richtingreconstructie van een diamantvormig station. zodat de driehoeken 134 en 123 met elkaar vergeleken
# kunnen worden.
# ook wordt op twee manieren de aankomsttijden berekent, een keer via 
# ProcessEvents._reconstruct_time_from_traces en via een weggeschreven tijden (prog: h5createtimestamps)
# er worden echter verschillen gevonden tussen beide methoden.
# deltatime_diamant bekijkt een station.

def func1(x, a, b, c): #fit met een gausskromme
    return a*np.exp(-((x-b)**2)/c**2)

def gauss(x, N, mu, sigma):
    return N*mlab.normpdf(x,mu,sigma)
    
def fitgauss(hulp, station, verschil):
    y, bins =  histogram(hulp, bins=np.arange(-41.25e-9, 41.25e-9, 2.5e-9), range=[-40e-9, 40e-9]) 
    x = (bins[:-1] + bins[1:])/2.
    N = max(y)
    
    f = lambda x, N, mu, sigma: N * scipy.stats.norm.pdf(x, mu, sigma)
    
    popt, pcov = curve_fit(gauss, x, y, [N * 1e-9, 0., 1e-9])
    xx = linspace(min(x), max(x), 1000)  
    versch = popt[1]
    waarden = 1e9*versch
    print 'van station %d en verschil %s' %(station, verschil) 
    print 'de verschuiving is %f'% waarden
    print 80*'-'
    # fig = plt.figure()    
    # ax = fig.add_subplot(111)
    # pylab.hist(hulp, bins=np.arange(-41.25e-9, 41.25e-9, 2.5e-9), range=[-40e-9, 40e-9], histtype='step', label='diff %d' % station)
     
    # pylab.plot(xx, gauss(xx, popt[0], popt[1], popt[2]), '--')
    # pylab.title('verschil in aankomsttijden voor %d \n' % station)
    # ax.annotate('tijdverschillen tussen %s van station %d' %(verschil, station), (-30e-9, 0.5*popt[0] ))
    # pylab.legend(loc='upper left' )
    # fig.show()
    
  
    
def deltatime_diamant(data, station):
    #events = data.getNode('/hisparc/cluster_amsterdam/station_%d/events' % station)
    
    events = data.getNode('/s%d/events' % station)# events = data.root.s502.events
    idx = events.getWhereList('(n1 > .5) & (n2 > .5) & (n3 > .5) & (n4 > .5)') 
 
    rows = events.readCoordinates(idx)
    
    timestamps = {'t1': [], 't2': [], 't3': [], 't4': []}
    
    t1 = 1.e-9 * events.readCoordinates(idx, field='t1')
    t2 = 1.e-9 * events.readCoordinates(idx, field='t2')
    t3 = 1.e-9 * events.readCoordinates(idx, field='t3')
    t4 = 1.e-9 * events.readCoordinates(idx, field='t4')
    
    dt1 = {'delta_2': []}
    dt12 = {'delta_12': []}
    dt13 = {'delta_13': []}
    dt14 = {'delta_14': []}
    
    dt1['delta_2'] = 0.333*(t1+t4+t3) - t2
    dt12['delta_12'] = t1 - t2
    dt13['delta_13'] = t1 - t3
    dt14['delta_14'] = t1 - t4
    
    hulp = dt12['delta_12']  
    verschil ='tijd 1-2'
    fitgauss(hulp, station, verschil)    
    
    hulp = dt13['delta_13']  
    verschil ='tijd 1-3'
    fitgauss(hulp, station, verschil) 

    hulp = dt14['delta_14']  
    verschil ='tijd 1-4'
    fitgauss(hulp, station, verschil)    
    


def angle_reconstruction(timestamp, cluster):
    # does not work, gives wrong phi and theta.
    
    station = cluster.stations[0]
       
    index, time = timestamp 
    c = 3.00e+8

    dt1 = time[0] - time[2]  
    dt2 = time[0] - time[3]  

    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)
    
    phi = arctan2((dt2 * 1e-9 * r1 * cos(phi1) - dt1 * 1e-9 * r2 * cos(phi2)),
                  (dt2 * 1e-9 * r1 * sin(phi1) - dt1 * 1e-9 * r2 * sin(phi2)) * -1)
    theta1 = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))
    theta2 = arcsin(c * dt2 * 1e-9 / (r2 * cos(phi - phi2)))

    # e1 = sqrt(self.rel_theta1_errorsq(theta1, phi, phi1, phi2, r1, r2))
    # e2 = sqrt(self.rel_theta2_errorsq(theta2, phi, phi1, phi2, r1, r2))

    # theta_wgt = (1 / e1 * theta1 + 1 / e2 * theta2) / (1 / e1 + 1 / e2)

    return theta1, theta2, phi

def deltatime(data):
    # parallellogram station, berekent tijdverschillen tussen zijden van het parallellogram 
    
    # events = data.root.hisparc.cluster_amsterdam.station_508.events
    events = data.root.s502.events
    idx = events.getWhereList('(n1 > .5) & (n2 > .5) & (n3 > .5) & (n4 > .5)') 
 
    rows = events.readCoordinates(idx)
    
    timestamps = {'t1': [], 't2': [], 't3': [], 't4': []}
    
    t1 = 1.e-9 * events.readCoordinates(idx, field='t1')
    t2 = 1.e-9 * events.readCoordinates(idx, field='t2')
    t3 = 1.e-9 * events.readCoordinates(idx, field='t3')
    t4 = 1.e-9 * events.readCoordinates(idx, field='t4')

           
    dt1 = {'dt12': [], 'dt43': []}
    dt2 = {'dt14': [], 'dt23': []}
    
    dt2['dt23'] = t2 - t3
    dt2['dt14'] = t1 - t4
    dt1['dt12'] = t1 - t2
    dt1['dt43'] = t4 - t3
        
    diff_time = {'deltat1': [], 'deltat2': []}
    for i in range(len(idx)):
        diff1 = []
        diff2 = []
        diff1 = dt1['dt12'][i]-dt1['dt43'][i]
        diff2 = dt2['dt14'][i]-dt2['dt23'][i]
        diff_time['deltat1'].append(diff1)
        diff_time['deltat2'].append(diff2)
    
    print diff_time['deltat1'][0:10]
    
    procestimegraph(diff_time)
    

def timeprocessEvents(timestamps, idx):
    
    
    dt1 = {'dt12': [], 'dt43': []}
    dt2 = {'dt14': [], 'dt23': []}
    
    for timestamp in timestamps:
        dt1['dt12'].append(timestamp[1][0]-timestamp[1][1])
        dt1['dt43'].append(timestamp[1][3]-timestamp[1][2])
        dt2['dt14'].append(timestamp[1][0]-timestamp[1][3])
        dt2['dt23'].append(timestamp[1][1]-timestamp[1][2])
    
    
    diff_time = {'deltat1': [], 'deltat2': []}
    for i in range(len(idx)):
        diff1 = []
        diff2 = []
        diff1 = dt1['dt12'][i]-dt1['dt43'][i]
        diff2 = dt2['dt14'][i]-dt2['dt23'][i]
        diff_time['deltat1'].append(diff1)
        diff_time['deltat2'].append(diff2)
    
    print diff_time['deltat1'][0:10]
    
    procestimegraph(diff_time)
    
    
def procestimegraph(diff_time):
    
    fig = plt.figure()    
    pylab.hist(diff_time['deltat1'], bins=np.arange(-41.25e-9, 41.25e-9, 2.5e-9), range=[-40e-9, 40e-9], histtype='step', label='diff 12-43')
    pylab.hist(diff_time['deltat2'], bins=np.arange(-41.25e-9, 41.25e-9, 2.5e-9), range=[-40e-9, 40e-9], 
    histtype='step', label='diff 14-23')
    
    hulp = np.array(diff_time['deltat1'])-np.array(diff_time['deltat2'])
    
    # pylab.hist(hulp, bins=np.arange(-41.25e-9, 41.25e-9, 2.5e-9), range=[-40e-9, 40e-9], histtype='step', label='verschil twee histogrammen')
    #pylab.yscale('log')
    pylab.title('parallellogram stations')
    pylab.legend(loc='upper left' )
    pylab.show()
    
    fig = plt.figure() 
    plt.axis('equal')
    plt.scatter(diff_time['deltat1'], diff_time['deltat2'])
    pylab.title('verschil in aankomsttijden')
    plt.xlim(-2e-7,2e-7)
    plt.ylim(-2e-7,2e-7)
    
    plt.show()
    
    
def find_row(data):
    events = data.root.s502.events
    
    idx = events.getWhereList('(n1 > .5) & (n2 > .5) & (n3 > .5) & (n4 > .5)')
    rows = events.readCoordinates(idx)
    
    return idx
    
def timings(data):
    #maakt tijden aan uit blobs
    idx = find_row(data)
    proces_time = ProcessEvents(data, '/s502')
    events = data.root.cluster_amsterdam.station_508.events
    rows = events.readCoordinates(idx)
    
    timings = [proces_time._reconstruct_time_from_traces(row) for row in rows]

    meetwaarden = zip(idx, timings)

    return meetwaarden, idx

if __name__ == '__main__':
    tekst = 'diamant10aug2013.h5'
    # data = tables.openFile('C:/Python/1tm4sept2013.h5', 'r')
    
    data = tables.openFile(os.path.join(DATAPATH, tekst), 'r')
    #deltatime(data)
    stations = [501, 502, 503, 504, 505, 506, 508, 509] 
    print 'metingen van verschil in tijd tussen detectoren van %s' %tekst
    
    for station in stations:
            
        deltatime_diamant(data, station) 
    
    
    # timestamps, idx = timings(data)   
    
    #timeprocessEvents(timestamps, idx)   
     
        
    rec_dir = DirectionReconstruction(data, min_n134=0.5)
    # rec_dir.station = data.root.core_reconstructions._v_attrs.cluster.stations[0]
    
    #row = data.root.s502.events
    # row = data.root.hisparc.cluster_amsterdam.station_508.events
    # idx = row.getWhereList('(n1 > .5) & (n2 > .5) & (n3 > .5) & (n4 > .5)')    
       
    # cluster = data.root.core_reconstructions._v_attrs.cluster
    # c_index = data.root.coincidences.c_index
    
    
    angles_134 = {'phi': [], 'theta': []}    
    angles_123 = {'phi': [], 'theta': []}
    
    
    
    #"""angle reconstructie in parallellogram station""""
    
    
    # for event in data.root.hisparc.cluster_amsterdam.station_508.events:
        # if (event['n1'] >= rec_dir.min_n134 and event['n2'] >= rec_dir.min_n134 and
            # event['n3'] >= rec_dir.min_n134 and event['n4'] >= rec_dir.min_n134):
            # detectors = [1,3,4]
            # theta, phi = rec_dir.reconstruct_angle(event, detectors)
            # angles_134['phi'].append(phi)
            # angles_134['theta'].append(theta)
            # detectors = [1,2,3]
            # theta, phi = rec_dir.reconstruct_angle(event, detectors)
            # angles_123['phi'].append(phi)
            # angles_123['theta'].append(theta)
        
    
    # fig = plt.figure()
    # pylab.hist(angles_123['phi'], bins=40, range=[-pi, pi], histtype='step', label='phi_123')
    # pylab.hist(angles_123['theta'], bins=40, range=[-pi, pi], histtype='step', label ='theta_123')
    # pylab.hist(angles_134['phi'], bins=40, range=[-pi, pi], histtype='step', label='phi_134')
    # pylab.hist(angles_134['theta'], bins=40, range=[-pi, pi], histtype='step', label='theta_134')
    # pylab.axis('auto')
   
    # pylab.xlabel(' Theta or Phi reconstructed')
    # pylab.ylabel('count')
    
    # pylab.title('Direction reconstruction (502)')
    # pylab.grid(True)
    # pylab.legend(loc='upper left' )
    # pylab.show()
    
    data.close()
    
    
    
    
    
