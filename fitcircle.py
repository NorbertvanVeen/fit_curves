import pylab as plt
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import fmin
import scipy
import tables
# Draw a fuzzy circle to test
n = 2 
j = 3
k = 5
X =[]
Y =[]
point4 = [83, 190.9]
point5 = [203.4,52.58]
point6 = [155.8, 211.17]
 
temp = []
coord =[]
with tables.openFile('C:/Python/station_coordinates.h5', 'r') as data2:
    detector = data2.root.coordinates.stations
    for detect in detector:
        temp = detect[5]
        coord.append(temp)
    data2.close()        

# X = (coord[n][0], coord[j][0], coord[k][0], point4[0], point5[0], point6[0])
# Y = (coord[n][1], coord[j][1], coord[k][1], point4[1], point5[1], point6[1])


X = [164.831, 207.14, 72.69, 171.721, 227.08]
Y = [209.758, 176.542, 178.124, 35.55, 107.044]
# N = 15
# THETA = np.random.random(15)*2*np.pi
# R     = 1.5 + (.1*np.random.random(15) - .05)
# X = R*np.cos(THETA) + 5
# Y = R*np.sin(THETA) - 2

# Choose the inital center of fit circle as the CM
xm = np.mean(X)
ym = np.mean(Y)

# Choose the inital radius as the average distance to the CM
cm = np.array([xm,ym]).reshape(1,2)
rm = cdist(cm, np.array([X,Y]).T).mean()

# Best fit a circle to these points
def err((w,v,r)):
    pts = [np.linalg.norm([x-w,y-v])-r for x,y in zip(X,Y)]
    return (np.array(pts)**2).sum()

xf,yf,rf = scipy.optimize.fmin(err,[xm,ym,rm])  

#ax = fig.add_subplot(1, 1, 1)
ax = plt.gca()
# Show the inital guess circle
#circ = plt.Circle((xm, ym), radius=rm, color='y',lw=2,alpha=.5)
#ax.add_patch(circ)

# Show the fit circle
circ = plt.Circle((xf, yf), radius=rf, color='b', alpha=0.1)
ax.add_patch(circ)
plt.scatter(X,Y, color ='g', zorder=4)
plt.show()
