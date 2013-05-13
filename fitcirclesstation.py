import numpy.linalg as linalg
import tables
import numpy as np
import matplotlib.pyplot as plt
import scipy 

from numpy import linspace
from scipy.spatial.distance import cdist
from scipy.optimize import fmin

from scipy import pi,sin,cos

def cost((x0, y0, a, b, phi)):
    Xe, Ye = f(x0, y0, a, b, phi)
    cost = 0
    distance_sq = (Xe - x0) ** 2 + (Ye - y0) ** 2
    cost += min(distance_sq)
    return cost

def optimize_parameters(center,axes, phi ):
    
    p0 = (center[0], center[1], axes[0], axes[1], phi)
    plt.plot(*f(*p0))

    xopt = fmin(cost, p0)
    print p0
    print xopt, cost(xopt)
    plt.plot(*f(*xopt))

    return xopt



def ellipse(axes,center, phi):
    Nb=50
    xpos,ypos=center
    radm,radn=axes
    an=phi

    co,si=cos(an),sin(an)
    the=linspace(0,2*pi,Nb)
    Xb=radm*cos(the)*co-si*radn*sin(the)+xpos
    Yb=radm*cos(the)*si+co*radn*sin(the)+ypos
    return Xb,Yb

def fitcircle(X,Y):
    xm = np.mean(X)
    ym = np.mean(Y)
    cm = np.array([xm,ym]).reshape(1,2)
    rm = cdist(cm, np.array([X,Y]).T).mean()

    xf,yf,rf = scipy.optimize.fmin(err,[xm,ym,rm])
    
    
    return xf, yf, rf
    
    
def err((w,v,r)):
    pts = [np.linalg.norm([x-w,y-v])-r for x,y in zip(X,Y)]
    return (np.array(pts)**2).sum()    
    
    
if __name__ == '__main__':
   
        
     
    
     # handmatige data volgorde van oplopende x is belangrijk! station =>503 504 506
    # X = [72.3569, 96.50, 156.634, 194.037, 201.198]
    # Y = [177.312, 201.243, 211.949, 45.6931, 184.535]
        
     # handmatige data volgorde van oplopende x is belangrijk! station =>501 502 505   
    X = [-249.032, -223.226, -172.903, -77.4194, 5.16]
    Y = [60.8856, -33.2103, -66.42, 171.587, 86.716]
    
    center, phi, axes = ellipse(X, Y)
    x1, y1 = ellipse(axes, center, phi)
    
    
    xf, yf, rf = fitcircle(X,Y)
    
    t = linspace(0, 2 * pi, 1000)
    f = lambda x0, y0, a, b, phi: \
        (x0 + a * cos(t) * cos(phi) - b * sin(t) * sin(phi),
         y0 + a * cos(t) * sin(phi) + b * sin(t) * cos(phi))

    
    xopt = optimize_parameters(center, axes, phi)
    
    print "fit van Norbert", center, axes, phi
    print "beste fit waarden voor stations =", xopt, stations
    
    circ = plt.Circle((xf, yf), radius=rf, color='k', alpha=0.1)
    ax = plt.gca()
    ax.add_patch(circ)
    plt.scatter(x1, y1, color ='k', zorder=4, label='punten op gefitte ellips')
    plt.scatter(X, Y, color ='g', zorder=4, label='punten op de boog')
    plt.legend()
    plt.show()
    