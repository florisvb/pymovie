import sys

sys.path.append('/home/floris/src/analysis')

from matplotlib.pyplot import figure, show
from matplotlib.patches import Ellipse, Rectangle
import numpy as np
import colorline
import matplotlib.pyplot as pyplot
import matplotlib.pyplot as plt
import scipy.linalg
import matplotlib.patches
import matplotlib.backends.backend_pdf as pdf
import pickle
import floris
import colorgrid
from matplotlib.patches import Patch
from matplotlib.patches import Arrow
from matplotlib.patches import Arc

import numpyimgproc as nim

import pygmovie as pm


def xy_kalman(npmovie):
    
    x = npmovie.kalmanobj.positions[:,1]
    y = npmovie.kalmanobj.positions[:,0]
    s = np.array([np.sqrt(npmovie.kalmanobj.positions[i,2]**2 + npmovie.kalmanobj.positions[i,3]**2) for i in range(len(npmovie.kalmanobj.positions))])
    
    cl = xy_trajectory(x,y,s)
    return cl
    
def xy_raw(npmovie, cl=None, colormap='jet'):

    x = npmovie.obj.positions[:,0]
    y = npmovie.obj.positions[:,1]
    s = np.array([np.sqrt(npmovie.obj.velocities[i,0]**2 + npmovie.obj.velocities[i,1]**2) for i in range(len(npmovie.obj.velocities))])

    if cl is None:
        cl = xy_trajectory(x,y,s, colormap=colormap)
        
    else:
        cl.colorline(x, y, s,colormap=colormap,linewidth=1)
    pyplot.show()
    return cl

def xy_trajectory(x, y, z, colorcode='s', norm=None, xlim=(0, 1024), ylim=(0,1024), figure=1, colormap='jet'):
    
    if norm is None:
        if colorcode == 's':
            norm = (0.02, .3)
    print colormap
    cl = colorline.Colorline(xlim=xlim, ylim =xlim, norm=norm, colormap = colormap, figure=figure)
    
    cl.colorline(x, y, z,linewidth=1)
    pyplot.show()
    
    return cl
    
    
    
    
def plot_movie_data(npmovie, movie):
    
    cl = xy_kalman(npmovie)
    
    # get strobe image:
    strobe_img = pm.strobe_image(movie)
    cl.ax0.imshow(strobe_img, pyplot.get_cmap('gray'))
    
    interval = 100.
    for i in range(len(npmovie.kalmanobj.positions)):
        if int(i/interval) == i/interval:
            center = npmovie.kalmanobj.positions[i]
            long_axis = npmovie.obj.long_axis[i]
            
            factor = 10
            dx = long_axis[1]*15
            dy = long_axis[0]*15
            arrow = Arrow(center[1], center[0], dx, dy, width=1.0, color='r')
            cl.ax0.add_artist(arrow)
            
    plt.show()
    return cl
        
        
        
    
    
    
