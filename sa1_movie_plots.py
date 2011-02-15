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
import classification as clas
import sa1_analysis as sa1

import copy

def get_active_frames(npmovie):
    # series of frames where position is not (0,0)
    try: 
        frame_of_landing = npmovie.frame_of_landing
    except:
        frame_of_landing = -1
        
    frames = np.nonzero(npmovie.obj.bool[0:frame_of_landing] > 100)[0]
    
    framediff = np.diff(frames)==1
    blobs = nim.find_biggest_blob(framediff)
    contframes = frames[ np.where(blobs==1)[0].tolist() ].tolist()

    return contframes
    
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
    calc_frame_of_landing(npmovie)
    frames = get_active_frames(npmovie)
    print frames
    time = np.array(frames)*1/float(npmovie.fps)
    cl = xy_kalman(npmovie, figure=figure, frames=frames)
    
    # get strobe image:
    strobe_img = strobe_from_npmovie(npmovie, interval=210, frames=frames)
    cl.ax0.imshow(strobe_img, pyplot.get_cmap('gray'))
    
    
    # axis parameters for subplots
    nxticks = 5
    nyticks = 3
    
    # subplot parameters
    n = 7
    h = (0.7-.05*(n-2))/float(n)
    
    subplots = []
    
    for s in range(n):
        ax = cl.fig.add_axes([0.71,0.1+(h+.05)*s,0.2,h])
        subplots.append(ax)
    
    xticks = np.linspace(0, 1, num=nxticks, endpoint=True).tolist()
    
    p = 0
    subplots[p].plot(time, npmovie.flycoord.worldangle[frames])
    subplots[p].set_xticks(xticks)
    default_yticks = subplots[p].get_yticks()
    yticks = np.linspace(default_yticks[0], default_yticks[-1], nyticks, endpoint=True).tolist()
    subplots[p].set_yticks(yticks)
    subplots[p].set_xlabel('time, seconds')
    subplots[p].set_ylabel('world angle')
    
    p += 1
    subplots[p].plot(time, npmovie.flycoord.postangle[frames])
    subplots[p].set_xticks(xticks)
    default_yticks = subplots[p].get_yticks()
    yticks = np.linspace(default_yticks[0], default_yticks[-1], nyticks, endpoint=True).tolist()
    subplots[p].set_yticks(yticks)
    subplots[p].set_xlabel('time, seconds')
    subplots[p].set_ylabel('angle to post')
    
    p += 1
    subplots[p].plot(time, npmovie.flycoord.slipangle[frames])
    subplots[p].set_xticks(xticks)
    default_yticks = subplots[p].get_yticks()
    yticks = np.linspace(default_yticks[0], default_yticks[-1], nyticks, endpoint=True).tolist()
    subplots[p].set_yticks(yticks)
    subplots[p].set_ylabel('slip angle')
    
    p += 1
    subplots[p].plot(time, npmovie.flycoord.velocities[frames,0])
    subplots[p].set_xticks(xticks)
    default_yticks = subplots[p].get_yticks()
    yticks = np.linspace(default_yticks[0], default_yticks[-1], nyticks, endpoint=True).tolist()
    subplots[p].set_yticks(yticks)
    subplots[p].set_ylabel('forward vel')
    
    p += 1
    subplots[p].plot(time, npmovie.flycoord.velocities[frames,1])
    subplots[p].set_xticks(xticks)
    default_yticks = subplots[p].get_yticks()
    yticks = np.linspace(default_yticks[0], default_yticks[-1], nyticks, endpoint=True).tolist()
    subplots[p].set_yticks(yticks)
    subplots[p].set_ylabel('sideways vel')
    
    p += 1
    subplots[p].plot(time, npmovie.kalmanobj.legs[frames])
    subplots[p].set_xticks(xticks)
    default_yticks = subplots[p].get_yticks()
    yticks = np.linspace(default_yticks[0], default_yticks[-1], nyticks, endpoint=True).tolist()
    subplots[p].set_yticks(yticks)
    subplots[p].set_ylabel('leg extension')
    
    p += 1
    subplots[p].plot(time, npmovie.flycoord.dist_to_post[frames])
    subplots[p].set_xticks(xticks)
    default_yticks = subplots[p].get_yticks()
    yticks = np.linspace(default_yticks[0], default_yticks[-1], nyticks, endpoint=True).tolist()
    subplots[p].set_yticks(yticks)
    subplots[p].set_ylabel('dist to post')
    

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



def pdf_movie_data(movie_dataset, scale = 10):
    
    # Initialize:
    pp =  pdf.PdfPages('sa1_movies_2_test.pdf')
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
        
        plt.Figure.set_figsize_inches(cl.fig, [2*scale,1*scale])
        plt.Figure.set_dpi(cl.fig, 72)
        
        cl.ax0.set_title(title)
        pp.savefig(f)
        plt.close(f)

    # Once you are done, remember to close the object:
    pp.close()
    
    
def reprocess_movies(movie_dataset, name=None, save=False):
    new_movie_dataset = {}
    for key,npmovie in movie_dataset.items():
        print 'processing: ', key
        
        pm.segment_fly(npmovie)
        pm.calc_obj_motion(npmovie)
        pm.smooth_legs(npmovie)
        clas.calc_fly_coordinates(npmovie)
        calc_frame_of_landing(npmovie)
        calc_timestamps(npmovie)
        
        try:
            align_sa1_flydra(npmovie)
        except:
            pass
            
        mnpmovie = pm.MiniNPM(npmovie)
        mnpmovie.behavior = copy.copy(npmovie.behavior)
        mnpmovie.path = copy.copy(npmovie.path)
        mnpmovie.posttype = copy.copy(npmovie.posttype)
        mnpmovie.extras = copy.copy(npmovie.extras)
        mnpmovie.id = copy.copy(npmovie.id)
        del(npmovie)
        
        new_movie_dataset.setdefault(mnpmovie.id, mnpmovie)
        
    if save:
        if name is None:
            name = 'sa1_movie_dataset_2_repro'
        pm.save(new_movie_dataset, name)
            
    return new_movie_dataset
        
def recalc_motion_movies(movie_dataset):
    for key,npmovie in movie_dataset.items():
        pm.calc_obj_motion(npmovie)
        pm.smooth_legs(npmovie)
        clas.calc_fly_coordinates(npmovie)
        calc_frame_of_landing(npmovie)
        
def calc_frame_of_landing(npmovie, threshold = 100.):
    frames = np.nonzero(npmovie.obj.bool > 20)[0].tolist()
    
    # if not landing behavior, return last frame as safety:
    if npmovie.behavior != 'landing':
        npmovie.frame_of_landing = frames[-1]
        print 'not landing'
        return npmovie.frame_of_landing
    
    # if landing, find frame where speed is below a threshold
    n_frames_landed = 0
    landing_frame = None
    for frame in frames[10:]:
        if npmovie.obj.wing_bool[frame] > threshold:
            n_frames_landed = 0
            landing_frame = None
        if npmovie.obj.wing_bool[frame] < threshold:
            if landing_frame is None:
                landing_frame = frame
            n_frames_landed += 1
        if n_frames_landed >= 50:
            npmovie.frame_of_landing = frame
            return npmovie.frame_of_landing
    
    
def calc_timestamps(npmovie):
    npmovie.timestamps = np.arange(0.,float(len(npmovie.uframes)), 1.)
    npmovie.timestamps *= 1./float(npmovie.fps)

def timestamp2frame(npmovie, timestamp):
    return float(npmovie.fps)*timestamp
        
#def align_sa1_flydra_fmin_func(dataset1, dataset2, index):
    
def interp_flydra_data_to_sa1(npmovie, data):
    flydra_raw_data = data
    flydra_raw_time = npmovie.trajec.fly_time
    flydra_interpolated_time = np.arange(flydra_raw_time[0], flydra_raw_time[-1], 1./float(npmovie.fps) )
    flydra_data = np.interp( flydra_interpolated_time, flydra_raw_time, flydra_raw_data )
    
    return flydra_data
    
def interp_sa1_data_to_flydra(npmovie, data, frames):
    sa1_raw_data = data
    sa1_raw_time = npmovie.timestamps[frames]
    sa1_interpolated_time = np.arange(sa1_raw_time[0], sa1_raw_time[-1], 1./float(100.) )
    sa1_data = np.interp( sa1_interpolated_time, sa1_raw_time, sa1_raw_data )
    
    return sa1_data
    
def convolve(a1, a2):
    
    if len(a1) > len(a2):
        along = a1
        ashort = a2
    else:
        along = a2
        ashort = a1
    
    conv = np.zeros( [len(along)-len(ashort)] )
    for i in range(len(conv)):    
        conv[i] = 1/np.sum(np.abs(along[i:i+len(ashort)] - ashort))
    
    return conv
    
    
def align_sa1_flydra(npmovie, plot=False):
    frames = get_active_frames(npmovie)
     
         
    sa1_data = np.diff(interp_sa1_data_to_flydra(npmovie, npmovie.flycoord.dist_to_post[frames].T[0], frames))
    flydra_data = np.diff(npmovie.trajec.dist_to_stim_r + npmovie.trajec.stimulus.radius)
    conv = convolve(sa1_data/sa1_data.max(), flydra_data/flydra_data.max())
    npmovie.sa1_start_index = np.argmax(conv)
            
    if plot:
        plt.plot(conv/np.max(conv))
        plt.plot(flydra_data/flydra_data.max())
        t = np.arange(npmovie.sa1_start_index, npmovie.sa1_start_index+len(sa1_data), 1)
        plt.plot(t, sa1_data/sa1_data.max()) 
            
        

    
    
    
    
