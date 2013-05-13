import datetime
import cPickle

import tables
import pylab
import numpy as np
from scipy.optimize import curve_fit
from scipy import ndimage
from pylab import histogram

MPV1, MPV2, MPV3, MPV4 = [], [], [], []

def fit_gauss(ph):
    for i in range(4):
		# bepaal de counts en de bingrootte in een range van 0 t/m 1000 ADC
        y, bins =  histogram(ph[i], bins=500, range=[0, 1000]) # loop for elke plaat
        x = (bins[:-1] + bins[1:])/2.  #let op neem bins om het midden te vinden van de bin
        
        if y[np.random.random_integers(200, 300)] != 0 : # ter controle of de data van de detector betrouwbaar is, aangezien plaat 4 503 geen goede data bevatte 3-01-2013
            # vanwege bartels (2012) kiezen we waar de ADC waarde tussen 150 en 410 zit voor de gauss fit, schat geeft de elementen waar tussen gefit moet worden.
            schat = np.where((x > 150) & (x < 600)) # geeft een tweeledig array terug    
            x1 = x[schat[0]] # geeft de waarde van de elementen gevonden met de determinatie hierboven in het eerste element
            #bovenstaande kan ook met x1 = fit_xa[np.where(x>150) & (x< 420)]
            y1 = y[schat] # bepaling van de count in de meting 
            max_min = ndimage.extrema(y1)
            max_y = ndimage.extrema(y1)[1] #maximale count y waarde piek van de gauss!)
        
            min_x = max_min[0]  #  de laagste waarde van de count
            max_x = max_min[1]  # hoogste waarde van de count -> a waarde in de gauss fit
            b_temp = max_min[3]
            b_place = b_temp[0] # b waarde
            b = x1[b_place] # de b-waarde van de gauss curve
        
            bound1 = max_x - (max_x - min_x)*0.75 
        
            if (max_x- min_x) <= 50:
                bound2 = max_x + (max_x - min_x)*1.5
            else:
                bound2 = max_x + (max_x - min_x)
        
            x2 = x1.compress((bound1 <= x1) & (x1 < bound2))
            y2 = y1.compress((bound1 <= x1) & (x1 < bound2))
        
              
            #popt, pcov = curve_fit(func1, x1, y1, [200, 250, 50] ) # werkende fit met een handmatige guess kleine verschillen. orde van 0.2 
            popt, pcov = curve_fit(func1, x1, y1, [max_x, b, 20] ) # fit met guess gebruikt werkt ook
        
            pylab.plot(func1(range(0,600), *popt), '--')
            peak = popt[1]
        
            if i == 0:
                MIP1 = peak
                MPV1.append(MIP1)  #bij herhalen  van de loop oor ander station append deze regel op de juiste manier?
            elif i == 1:
                MIP2 = peak
                MPV2.append(MIP2)
            elif i == 2:
                MIP3 = peak
                MPV3.append(MIP3)
            elif i == 3:
                MIP4 = peak
                MPV4.append(MIP4)
        
            print 'The MPV of detector',i + 1, 'lies at', peak, 'ADC'
        
        else:
            
            print 'The data of the detector ',i + 1, 'could not be fitted to a gauss curve'
    
    
def func(x, a, b, c): #fit met een parabool
    return a * (x - b) ** 2 + c   
    
def func1(x, a, b, c): #fit met een gausskromme
    return a*np.exp(-((x-b)**2)/c)
    
def hist_ph(ph):
    for i in range(4):
        pylab.hist(ph[i], bins=500, range=[0, 1000], histtype='step')
        pylab.axis('auto')
        pylab.xlabel('pulseheight (ADC)')
        pylab.ylabel('count')
        pylab.yscale('log')
        pylab.title('Pulseheight histogram (501)')
        pylab.grid(True)
        pylab.show()


if __name__ == '__main__':

    #get_coincidences()
    # in deze file staat nog maar een station, meerdere hebben!
    # pkl_file = open('datatest_coin.pkl', 'rb')
    # pulseheights = cPickle.load(pkl_file)
    data = tables.openFile('C:/Python/Datasets/datatest1.h5', 'a')
    pulseheights = data.root.s503.events.col('pulseheights')
    pulseheights = np.array(pulseheights)
    pulseheights = zip(*pulseheights)
    hist_ph(pulseheights)
      
    # indexltr =['A', 'B', 'C', 'D']
    # for i in range(len(indexltr)):
        # pulserij = pulseheights[i]
        # print indexltr[i]
     
    fit_gauss(pulseheights)
    
    # plot_ph_hist(pulseheights)