from numpy import nan, pi, sin, cos, tan, arctan2, arctan, arcsin
import numpy as np
import tables
from math import *
import scipy
import pylab
import matplotlib.pyplot as plt
import operator


from scipy.stats import poisson, chi2
from scipy.special import gamma
from pylab import histogram
from scipy.optimize import curve_fit
from itertools import product, combinations
from scipy.optimize import leastsq

from fit_form import coord_tcc, coin_stations


def hist_Ne(sel_coinc_ids, coords, data, n):
    global histNe2
    c_index = data.root.coincidences.c_index
    observ = data.root.coincidences.observables
    core_rec = data.root.core_reconstructions.reconstructions
    
    histNe = core_rec.readCoordinates(coords, field='reconstructed_shower_size')
    #histNe = [x for x in histNe if x > 0]  # for showersize smaller than 0 
    
    d = 10**10.2
    
    histNe2 = [x*d for x in histNe]
    #histNe *= 10 ** 10.2
    
      
    pylab.hist(np.log10(histNe2), 100, log=True) # histtype="step"
       
    pylab.xlabel('Showerenergy log(eV)')
    pylab.ylabel('count')
    pylab.title('showersize bij N==%s' %n)
    pylab.ylim(ymin=1)
    pylab.grid(True)
    pylab.show()
    
    return histNe2
    
def func(x, a, b):
    return a*(10**bx)
    
def calc_slope(histN):
    y, x = np.histogram(np.log10(histN), 100)  #y = count x = bins
    x = x[:-1]
    index, value = max(enumerate(y), key=operator.itemgetter(1))
    x1 = x[index:]
    y1 = y[index:]
    x0 = x[index]
    y0 = y[index]
    y2 = y[index+30]
    x2 = x[index+30]
    
    print x1
    popt, _ = scipy.optimize.curve_fit(lambda x, a, b: a*(10**x)**b, x1, y1, [400,-5])
    a,bc = popt
   
            
    slope,intercept=np.polyfit(y1, x1,1)
    #b = y0 - (10**slope)*x0
    b = (np.log(y0)-np.log(y2))/(np.log(x0)-np.log(x2))
    #or calculate:
    yz =func(x1, a, bc)
    print yz
    plt.scatter(x1, yz,color ='red', zorder=4, label='linear fit power %s'%a)
    # test = [400. * (10.**x)**-3 for x in range(15, 20, 0.2)]
    # plt.plot(np.arange(15, 20, 0.2), test)
    plt.legend()
    plt.show()
    
    print "helling", slope 
    print "macht", b
    print "helling tweede manier", a, bc
 

    
if __name__ == '__main__':
    data = tables.openFile('C:/Python/results-feb-2012.h5', 'r')
    
    n = 4
    tcc = 10
    coinc_ids, tccs, densities = coin_stations(data,n)
    sel_coinc_ids, result_coinc_ids, coords = coord_tcc(coinc_ids, tccs, data,tcc)
    histN = hist_Ne(sel_coinc_ids, coords, data, n)
    calc_slope(histN)