#!/usr/bin/env python

import mdtraj as md
import numpy as np
import scipy.interpolate
import math as math
from scipy.ndimage import gaussian_filter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
from pylab import *
#import seaborn as sns
from scipy.misc import derivative

from matplotlib import pyplot as plt
from matplotlib import animation

fname="fes_2d"

# first read the PMF from fes
# find reference point as fes = 0

ref_x = [0.6, 2.4]
ref_y = [0.9, 2.7]

def readpmf (name, x1, y1, z1, xr, yr, maxval=100):
    x, y, z = np.genfromtxt(name, unpack=True, comments="#", usecols=(0,1,2))
    minz=3000.0
    minx=3000.0
    miny=3000.0
    for i in range(len(z)):
        if ((np.isinf(z[i]) or z[i]>maxval)):
            z[i]=maxval

        x1=np.append(x1,x[i])
        y1=np.append(y1,y[i])
        z1=np.append(z1,z[i])

        if (x[i] > ref_x[0] and x[i] < ref_x[1] and y[i] > ref_y[0] and y[i] < ref_y[1]) :
            if (z[i] < minz) :
                minz=z[i]
                minx=x[i]
                miny=y[i]

    z1-=minz
    print("The min is %f" % minz)
    print("The min(x,y) is %f, %f" % (minx, miny))
    f_in = open(name,'r')
    f_in.readline()
    #f_in.readline()
    xr[0] = float(f_in.readline().split()[3])
    xr[1] = float(f_in.readline().split()[3])
    f_in.readline()
    f_in.readline()
    yr[0] = float(f_in.readline().split()[3])
    yr[1] = float(f_in.readline().split()[3])
    f_in.close()
    print(xr,yr)
    #z1=np.divide(z1,4.184) # convert to kcal/mol from kJ/mol
    x1=np.multiply(x1,1.0)
    y1=np.multiply(y1,1.0)
    return x1,y1,z1,minz,xr,yr

xa=np.array([])
ya=np.array([])
za=np.array([])

maxval=100 #kJ/mol
nbin=191 #number of bins

BASE="fes_2d"
FES_FILE=".dat"
xa,ya,za,minfe,xr,yr=readpmf(BASE+FES_FILE,xa,ya,za,maxval=maxval,xr=np.array([0,104]),yr=np.array([0,104]))
xi1,yi1=np.linspace(xa.min(),xa.max(),nbin),np.linspace(ya.min(),ya.max(),nbin)
xi1 = xi1

xx,yy = np.meshgrid(xi1, yi1)

# this is where I would average fes's but I'm just using 1 fes.dat for this example
z1 = za
z1sq = np.square(za)
x1 = xa
y1 = ya

zi1 = scipy.interpolate.griddata((x1,y1),z1,(xx,yy),method='cubic')
zi1 = gaussian_filter(zi1, sigma=1.5) 
zi1sq = scipy.interpolate.griddata((x1,y1),z1sq,(xx,yy),method='cubic')
zi1sq = gaussian_filter(zi1sq, sigma=1.5)

for i in range(len(zi1)):
  for j in range(len(zi1[0])):
    if (math.isnan(zi1[i][j])):
       zi1[i][j]=maxval
    if (math.isnan(zi1sq[i][j])):
       zi1sq[i][j]=maxval

ziter=scipy.interpolate.interp2d(xi1,yi1,zi1,kind='cubic')
zsqiter=scipy.interpolate.interp2d(xi1,yi1,zi1sq,kind='cubic')
zdrvxdata=ziter(xi1,yi1,1,0)
zdrvydata=ziter(xi1,yi1,0,1)

#parameterize the partial derivative function f(x,y)=dz/dx
zdrvx=scipy.interpolate.interp2d(xi1,yi1,zdrvxdata,kind='cubic')
zdrvy=scipy.interpolate.interp2d(xi1,yi1,zdrvydata,kind='cubic')


# INITIALIZE STRING METHOD HERE --------------------------------------------
# from point a to e in this case
xa= 0.54
ya= 0.83

xb= 1.01
yb= 1.28

xc= 1.48
yc= 1.735

xd= 1.95
yd= 2.19

xe= 2.42
ye= 2.64

#number of points 
n1=40
n1sub=10 #4 lines, num points per line

#splice the intial path to n1 equidistant points
g1=np.linspace(0,1,n1)
g1sub=np.linspace(0,1,n1sub)
g2sub=np.linspace(0,1,n1sub+1)

x1v = (xb-xa)*g1sub+xa
y1v = (x1v-xa)*(yb-ya)/(xb-xa)+ya

x2v = (xc-xb)*g2sub+xb
y2v = (x2v-xb)*(yc-yb)/(xc-xb)+yb

x3v = (xd-xc)*g2sub+xc
y3v = (x3v-xc)*(yd-yc)/(xd-xc)+yc

x4v = (xe-xd)*g2sub+xd
y4v = (x4v-xd)*(ye-yd)/(xe-xd)+yd

x = np.concatenate((x1v,x2v[1:],x3v[1:],x4v[1:]), axis=0)
y = np.concatenate((y1v,y2v[1:],y3v[1:],y4v[1:]), axis=0)

dx = x - np.roll(x,1)
#normalize the scale of y distances
dy= (y - np.roll(y,1))/(yi1.max()-yi1.min())*(xi1.max()-xi1.min())

dx[0]=0
dy[0]=0

#length of the path
sq=np.sqrt(np.square(dx)+np.square(dy))

lxy = np.cumsum(sq)
lxy = lxy/lxy[n1-1];

#splice the intial path to equidistant points
x=np.interp(g1,lxy,x)
y=np.interp(g1,lxy,y)


# STRING CONVERGENCE METRICS ---------------------------------------

maxstep=20000
#maxstep=0
hx = 0.0000005 #may need to adjust this by orders of magnitude (e.g. maybe 1e-7 necessary)
hy = 0.0000005 #each variable can also have different step sizes they need
#tolerance
tol1 = 6e-7 #also adjust based on how big the gradients are
tol = 0.1

for k in range(maxstep):

    # calculation of the x and y-components of the force, dVx and dVy respectively
    print("Iteration %d, err %9.8f" %(k,tol), end='\r')
    dVx=np.array([])
    dVy=np.array([])
    for j in range(len(x)):
        dVx=np.append(dVx,zdrvx(x[j],y[j])[0])
        dVy=np.append(dVy,zdrvy(x[j],y[j])[0])
    x0 = x
    y0 = y

    # string steps:

    # 1. evolve
    x = x - hx*dVx
    y = y - hy*dVy

    # 2. reparametrize
    dx = x - np.roll(x,1)
    dy = (y - np.roll(y,1))/(yi1.max()-yi1.min())*(xi1.max()-xi1.min())
    dx[0]=0
    dy[0]=0
    sq=np.sqrt(np.square(dx)+np.square(dy))

    lxy = np.cumsum(sq)
    lxy = lxy/lxy[n1-1];
    x=np.interp(g1,lxy,x)
    y=np.interp(g1,lxy,y)

    tol = (np.linalg.norm(x-x0)+np.linalg.norm(y-y0))/n1

    if (tol <= tol1):
        break
    #print(tol)

if (tol > tol1):
  print('The calculation failed to converge after {0:5d} iterations\n'.format(maxstep))
else:
  print('The calculation converged after {0:5d} iterations\n'.format(k+1))

dx = x - np.roll(x,1)
dy = (y - np.roll(y,1))/(yi1.max()-yi1.min())*(xi1.max()-xi1.min())
dx[0]=0
dy[0]=0
sq=[np.sqrt(dx[i]**2+dy[i]**2) for i in range(len(dx))]
lxy = np.cumsum(sq)
lxy = lxy/lxy[n1-1];

# 3. higher resolution interpolation
n2=200
g2=np.linspace(0,1,n2)

x=np.interp(g2,lxy,x)
y=np.interp(g2,lxy,y)

font = { 'size'   : 18}
axes = {'linewidth': 2,
        'labelsize': 'large',
       }
#xtick = {'major.size' : 10,
#         'minor.size' : 5,
#         'major.width' : 2,
#         'minor.width' : 2
#        }
#ytick=xtick
matplotlib.rc('font',**font)
#matplotlib.rc('axes',**axes)
#matplotlib.rc('xtick',**xtick)
#matplotlib.rc('ytick',**ytick)
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'


def plotfes(name,xi,yi,zi,z1,x,y):
    fig=plt.figure(figsize=(11,8),dpi=600)
    ax=fig.add_axes([0.25,0.18,0.65,0.70])

    ax.plot(x,y,color='r',lw=2,linestyle='--')
    ax.plot(x,y,color='r',marker='.',markersize=16,linestyle='None',markevery=10)
    ax.set_xlabel(r"$Distance_{0}$ (nm)",fontsize=25, labelpad=10)
    ax.set_ylabel(r"$Distance_{1}$ (nm)",fontsize=25, labelpad=10)

    majorLocator=MultipleLocator(10)
    majorFormator=FormatStrFormatter('%2d')
    minorLocator=MultipleLocator(2)

    #ax.xaxis.set_major_locator(majorLocator)
    #ax.xaxis.set_major_formatter(majorFormator)
    #ax.xaxis.set_minor_locator(minorLocator)
    #ax.tick_params(axis='x', pad=8)

    #majorLocator=MultipleLocator(10.0)
    #majorFormator=FormatStrFormatter('%0.2f')
    #minorLocator=MultipleLocator(2.0)
    #ax.yaxis.set_major_locator(majorLocator)
    #ax.yaxis.set_major_formatter(majorFormator)
    #ax.yaxis.set_minor_locator(minorLocator)
    #ax.tick_params(axis='y', pad=10)
  
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(18)
    ax.set_xlim(0.55, 2.42)
    ax.set_ylim(0.85, 2.72)
    
    ax.set_xticks(np.arange(0.60,2.451,0.4))
    ax.set_yticks(np.arange(0.90,2.751,0.4))
    
    #ax.set_ylim(0.2, 1.0)
    cmap=plt.get_cmap('bone')

    vcap = 1.1
    im=ax.imshow(zi,interpolation='bilinear',origin='lower',cmap=cmap,extent=[x1.min(),x1.max(),y1.min(),y1.max()],aspect='auto') #,vmin=z1.min(),vmax=vcap
    #levels=[x for x in np.arange(-10.0,vcap,0.1)]  
    #levels2=[x for x in np.arange(-10.0,vcap,1.0)]
    #im2=ax.contourf(xi,yi,zi,levels=levels,cmap=cmap)

    #im3=ax.contour(xi,yi,zi,levels=levels2,colors='k',linewidth=4)

    #add colorbar label
    cb2=fig.colorbar(im,ax=ax,format='%d',label='',fraction=0.185,pad=0.08)
    cb2.set_label("Free Energy [kcal/mol]", fontsize=30, rotation=270, labelpad=40)
    #if (name == 'diff.pdf') :
    #    cb2.set_ticks([int(x) for x in np.arange(z1.min(),z1.max(),2)])
    #else:
    #    cb2.set_ticks([int((x)/1)*1 for x in np.arange(0,vcap,2)])
    #cb2.ax.tick_params(labelsize=30)
    fig.savefig(name)


plotfes('path_fes.png',xx,yy,zi1,z1,x,y)

z1d=np.array([])
dis_xy_tol=np.sqrt((x[len(x)-1]-x[0])**2+(y[len(y)-1]-y[0])**2)

output=open('min_path.dat',"w")
for i in range(len(x)):
    # Add this line to compute the distance on the line
    dis_xy=np.sqrt((x[i]-x[0])**2+(y[i]-y[0])**2)
    dis_xy=dis_xy/dis_xy_tol
    E2 = zsqiter(x[i],y[i])[0]
    E1 = ziter(x[i],y[i])[0]
    stdev = ( abs(E2 - E1*E1) )**0.5
    stderr = stdev / 1.0 #only 1 replica in this example
    output.write('{0:14.8f}  {1:14.8f}  {2:14.8f}  {3:14.8f} {4:14.8f}\n'.format(x[i],y[i],dis_xy,ziter(x[i],y[i])[0],stderr))
