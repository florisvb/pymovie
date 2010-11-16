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
import pickle
import pygmovie as pm

import copy

def load(filename):
    fname = (filename)  
    fd = open( fname, mode='r' )
    npmovie = pickle.load(fd)
    fd.close()
    return npmovie
    
    
def xy_kalman(npmovie):
    
    x = npmovie.kalmanobj.positions[:,1]
    y = npmovie.kalmanobj.positions[:,0]
    s = npmovie.kalmanobj.speed
    
    cl = xy_trajectory(x,y,s)
    return cl
    
def xy_raw(npmovie, cl=None, colormap='jet', strobe=False, movie=None):

    x = npmovie.obj.positions[:,1]
    y = npmovie.obj.positions[:,0]
    s = np.array([np.sqrt(npmovie.obj.velocities[i,0]**2 + npmovie.obj.velocities[i,1]**2) for i in range(len(npmovie.obj.velocities))])

    if cl is None:
        cl = xy_trajectory(x,y,s, colormap=colormap)
        
    else:
        cl.colorline(x, y, s,colormap=colormap,linewidth=1)
    #pyplot.show()

    if strobe is True:
        strobe_img = pm.strobe_image(movie)
        cl.ax0.imshow(strobe_img, pyplot.get_cmap('gray'))

    return cl

def xy_trajectory(x, y, z, colorcode='s', norm=None, xlim=(0, 1024), ylim=(0,1024), figure=1, colormap='jet'):
    
    if norm is None:
        if colorcode == 's':
            norm = (0.02, .3)
    print colormap
    cl = colorline.Colorline(xlim=xlim, ylim =xlim, norm=norm, colormap = colormap, figure=figure)
    
    cl.colorline(x, y, z,linewidth=1)
    #pyplot.show()
    
    return cl
    
    
    
    
def plot_movie_data(npmovie):
    
    cl = xy_kalman(npmovie)
    
    # get strobe image:
    strobe_img = strobe_from_npmovie(npmovie)
    cl.ax0.imshow(strobe_img, pyplot.get_cmap('gray'))
    
    interval = 100.
    for i in range(len(npmovie.kalmanobj.positions)):
        if int(i/interval) == i/interval:
            center = npmovie.kalmanobj.positions[i]
            long_axis = npmovie.kalmanobj.long_axis[i]
                       
            factor =  npmovie.obj.axis_ratio[i][0]*4.
            factor = min(15., factor)
            print i, factor
            
            dx = long_axis[1]*factor
            dy = long_axis[0]*factor
            arrow = Arrow(center[1], center[0], dx, dy, width=1.0, color='r')
            cl.ax0.add_artist(arrow)
            
    #plt.show()
    return cl
    
    
    
def plot_movie_dict(movie_dict, cl=None):
    
    x = movie_dict['KalmanObj'].positions[:,1]
    y = movie_dict['KalmanObj'].positions[:,0]
    s = movie_dict['KalmanObj'].speed
    cl = xy_trajectory(x,y,s)
    
    # get strobe image:
    strobe_img = movie_dict['StrobeImg']
    cl.ax0.imshow(strobe_img, pyplot.get_cmap('gray'))
    
    interval = 100.
    for i in range(len(movie_dict['KalmanObj'].positions)):
        if int(i/interval) == i/interval:
            center = movie_dict['KalmanObj'].positions[i]
            long_axis = movie_dict['KalmanObj'].long_axis[i]
                       
            factor =  movie_dict['Obj'].axis_ratio[i][0]*4.
            factor = min(15., factor)
            print i, factor
            
            dx = long_axis[1]*factor
            dy = long_axis[0]*factor
            arrow = Arrow(center[1], center[0], dx, dy, width=1.0, color='r')
            cl.ax0.add_artist(arrow)
            
    #plt.show()
    return cl
    
def strobe_from_npmovie(npmovie):

    bkgrd = npmovie.background
    strobe_img = copy.copy(bkgrd)
    
    interval = 10
    i = 0
    while i < len(npmovie.uframes):
        if npmovie.uframes[i].uimg is not None:
            
            uimg = npmovie.uframes[i].uimg
            indices = npmovie.uframes[i].indices
            center = npmovie.uframes[i].center
            blank = 255*np.ones_like(bkgrd)
            #blank.dtype = float
            
            blank[indices[0]:indices[1], indices[2]:indices[3]] = uimg
            strobe_img = nim.darken(strobe_img, blank)
                
        i += interval
        
    return strobe_img
    
def plot_trajectory(movie_info, movie, strobe=True, colorcode='s', norm=None, xlim=(0, 1024), ylim=(0,1024), figure=None, colormap='jet'):
    
    cl = colorline.Colorline(xlim=xlim, ylim =xlim, norm=norm, colormap = colormap, figure=figure)
    
    movie_dict = movie_info[movie]
    
        
    x = movie_dict['KalmanObj'].positions[:,1]
    y = movie_dict['KalmanObj'].positions[:,0]
    s = movie_dict['KalmanObj'].speed
    cl.colorline(x,y,s, linewidth=1)
    
    # get strobe image:
    strobe_img = movie_dict['StrobeImg']
    cl.ax0.imshow(strobe_img, pyplot.get_cmap('gray'))
    
    interval = 100.
    for i in range(len(movie_dict['KalmanObj'].positions)):
        if int(i/interval) == i/interval:
            center = movie_dict['KalmanObj'].positions[i]
            long_axis = movie_dict['KalmanObj'].long_axis[i]
                       
            factor =  movie_dict['Obj'].axis_ratio[i][0]*4.
            factor = min(15., factor)
            #print i, factor
            
            dx = long_axis[1]*factor
            dy = long_axis[0]*factor
            arrow = Arrow(center[1], center[0], dx, dy, width=1.0, color='r')
            cl.ax0.add_artist(arrow)
            
    title = movie
    cl.ax0.set_title(title)
    cl.ax0.set_xlabel('pixels')
    cl.ax0.set_ylabel('pixels')
                    
                    
            
def pdf_landings(movie_info):
    
    plt.close('all')
    
    pp =  pdf.PdfPages('SA1_landings.pdf')
    f = 0
    
    for movie in movie_info.keys():
        if movie_info[movie]['Behavior'] == 'landing':
            # As many times as you like, create a figure fig, then either:
            f += 1
            plot_trajectory(movie_info, movie, figure = f)
            pp.savefig(f)
            plt.close(f)
        
    
    
    # Once you are done, remember to close the object:
    pp.close()
    print 'closed'
            
    
    
        
        
        
    
    
    
