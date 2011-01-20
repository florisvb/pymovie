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

def get_active_frames(npmovie):
    # series of frames where position is not (0,0)

    p = 0
    f = 0
    interval = 10
    while p < 1:
        f += 1
        tmp = [np.sum(npmovie.kalmanobj.positions[ff]) for ff in range(f,f+interval)]
        p = np.min( tmp )
    start = f
    
    p = 2
    f = start
    interval = 2
    while p > 1:
        f += 1
        if f+interval > len(npmovie.uframes)-1:
            break
        tmp = [np.sum(npmovie.kalmanobj.positions[ff]) for ff in range(f,f+interval)]
        p = np.min( tmp )
    stop = f
        
    active_range = np.arange(start, stop).tolist()
    return active_range
    
def get_all_frames(npmovie):
    
    active_range = np.nonzero(npmovie.kalmanobj.positions[:,0])[0].tolist()
    
    return active_range
        
def load(filename):
    fname = (filename)  
    fd = open( fname, mode='r' )
    npmovie = pickle.load(fd)
    fd.close()
    return npmovie
    
    
def xy_kalman(npmovie, frames=None, figure=None):
    
    if frames is None:
        x = npmovie.kalmanobj.positions[:,1]
        y = npmovie.kalmanobj.positions[:,0]
        s = npmovie.kalmanobj.speed[:,0]

    else:
        x = npmovie.kalmanobj.positions[frames,1]
        y = npmovie.kalmanobj.positions[frames,0]
        s = npmovie.kalmanobj.speed[frames,0]
        
    cl = xy_trajectory(x,y,s,figure=figure)
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

def xy_trajectory(x, y, z, colorcode='s', norm=None, xlim=(0, 1024), ylim=(0,1024), figure=None, colormap='jet'):
    
    if norm is None:
        if colorcode == 's':
            norm = (0.02, .3)
    print colormap
    cl = colorline.Colorline(xlim=xlim, ylim =xlim, norm=norm, colormap = colormap, figure=figure, hide_colorbar=True)
    
    cl.colorline(x, y, z,linewidth=1)
    #pyplot.show()
    
    return cl
    
    
    
    
def plot_movie_data(npmovie, show_wings=False, figure=None, legthresh=50):


    frames = get_all_frames(npmovie)
    cl = xy_kalman(npmovie, figure=figure, frames=frames)
    
    # get strobe image:
    strobe_img = strobe_from_npmovie(npmovie, interval=210, frames=frames)
    cl.ax0.imshow(strobe_img, pyplot.get_cmap('gray'))
    
    interval = 70
    i = frames[0]
    while i < frames[-1]:
    
        # plot body orientation vector
        center = npmovie.kalmanobj.positions[i]
        long_axis = npmovie.kalmanobj.long_axis[i]
                   
        factor =  npmovie.obj.axis_ratio[i][0]*4.
        factor = min(15., factor)
        
        dx = long_axis[1]*factor
        dy = long_axis[0]*factor
        
        legs = npmovie.kalmanobj.legs[i]
        if legs > legthresh:
            arrow = Arrow(center[1], center[0], dx, dy, width=1.0, color='r')
        else:
            arrow = Arrow(center[1], center[0], dx, dy, width=1.0, color='b')
    
        cl.ax0.add_artist(arrow)
        
        if show_wings:
            # plot wing orientation vectors
            wingR = npmovie.kalmanobj.wingcenterR[i]
            if not np.isnan(wingR[0]):
                arrow = Arrow(center[1], center[0], wingR[1]-30, wingR[0]-30, width=1.0, color='b')
                cl.ax0.add_artist(arrow)
            
            wingL = npmovie.kalmanobj.wingcenterL[i]
            if not np.isnan(wingL[0]): 
                arrow = Arrow(center[1], center[0], wingL[1]-30, wingL[0]-30, width=1.0, color='b')
                cl.ax0.add_artist(arrow)
        
        
        i += interval
            
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
    
def strobe_from_npmovie(npmovie, interval=200, frames=None):

    bkgrd = npmovie.background
    strobe_img = copy.copy(bkgrd)
    
    i = 0
    
    if frames is None:
        frames = np.arange(0,len(npmovie.uframes),interval)
    else:
        frames = np.arange(frames[0], frames[-1], interval)    
    
    for i in frames:
        if npmovie.uframes[i].uimg is not None:
            
            uimg = npmovie.uframes[i].uimg
            indices = npmovie.uframes[i].indices
            center = npmovie.uframes[i].center
            blank = 255*np.ones_like(bkgrd)
            #blank.dtype = float
            
            blank[indices[0]:indices[1], indices[2]:indices[3]] = uimg
            strobe_img = nim.darken(strobe_img, blank)
                
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
            
    
    
        
        
        
    
############## wing kinematics  #####################

def plot_wingbeat_frequency(npmovie):
    frames = get_active_frames(npmovie)
    
    waR = np.reshape(npmovie.kalmanobj.wingangleR[frames], (len(frames)))
    interval = 1000

    i = 0
    seg = waR[i:i+interval]
    sp = np.fft.fft(seg)
    freq = np.fft.fftfreq(len(seg))
    plt.plot(freq*npmovie.fps, np.abs(sp.real))
    
def plot_wingbeats(npmovie):

    # REQUIRES that you have run sa1a.get_flydra_frame()
    
    frames = get_active_frames(npmovie)
    t = np.array(frames)*1./npmovie.fps
    
    waR = np.reshape(npmovie.kalmanobj.wingangleR[frames], (len(frames)))
    waL = np.reshape(npmovie.kalmanobj.wingangleL[frames], (len(frames)))
     
    plt.subplot(421)
    plt.plot(t, waR)
    plt.xlabel('time')
    plt.ylabel('right wing angle (radians)')

    plt.subplot(423)
    plt.plot(t, waL)
    plt.xlabel('time')
    plt.ylabel('left wing angle (radians)')
    
    r = npmovie.kalmanobj.polarpos[frames][:,0]
    theta = npmovie.kalmanobj.polarpos[frames][:,1]
    
    
    z = np.zeros_like(r)
    i = 0
    for uframe in npmovie.uframes[frames[0]:frames[-1]]:
        tmp = npmovie.trajec.positions[uframe.flydraframe,2] 
        z[i] = tmp
        i += 1
            
    plt.subplot(422)
    plt.plot(t, r)
    plt.xlabel('time')
    plt.ylabel('fly position, dist to post, pixels')
    
    plt.subplot(424)
    plt.plot(t, theta)
    plt.xlabel('time')
    plt.ylabel('angle to post, radians')
    
    if 1:
        plt.subplot(426)
        plt.plot(t, z)
        plt.xlabel('time')
        plt.ylabel('elevation, meters')
        
    plt.subplot(425)
    plt.plot(t, npmovie.kalmanobj.legs[frames])
    plt.xlabel('time')
    plt.ylabel('leg extension, arb. units')
    
    plt.subplot(427)
    plt.plot(t, npmovie.kalmanobj.axis_ratio[frames])
    plt.xlabel('time')
    plt.ylabel('body length ratio (pitch est.)')
    
    plt.subplot(428)
    plt.plot(t, npmovie.kalmanobj.speed[frames])
    plt.xlabel('time')
    plt.ylabel('speed, pixels/second')
    
    plt.show()



def pdf_movie_data(movie_dataset):
    
    # Initialize:
    pp =  pdf.PdfPages('sa1_movies_2.pdf')
    f = 0
    
    for key,mnpmovie in movie_dataset.items():
        f += 1
        cl = plot_movie_data(mnpmovie, show_wings=False, figure=f)
        
        try:
            extras = mnpmovie.extras[0]
        except:
            extras = '' 
        print mnpmovie.id, mnpmovie.behavior, extras
        if extras == 'none':
            extras = ''
        title = mnpmovie.id + ' ' + mnpmovie.behavior + ' ' + extras
        
        plt.Figure.set_figsize_inches(cl.fig, [20,20])
        plt.Figure.set_dpi(cl.fig, 300)
        
        cl.ax0.set_title(title)
        pp.savefig(f)
        plt.close(f)

    # Once you are done, remember to close the object:
    pp.close()
    
    
    

    
