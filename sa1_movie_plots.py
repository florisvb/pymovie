import sys
sys.path.append('/home/floris/src/analysis')
sys.path.append('/home/floris/src/floris')
sys.path.append('/home/floris/src/fsee_utils/src/floris')

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

import adskalman.adskalman as adskalman
import numpyimgproc as nim
import pickle
import pygmovie as pm
import classification as clas
import sa1_analysis as sa1

import kmeans_npmovie as kn

import scipy.misc

import camera_math
import fsee_timeseries

import flydra.reconstruct

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
    
    
def xy_kalman(npmovie, frames=None, figure=None, colormap='jet'):
    
    if frames is None:
        frames = get_all_frames(npmovie)

    posrot = sa1_to_flydra_data_transform(npmovie.kalmanobj.positions[frames,:], fix_sign=True)
    x = posrot[:,0]
    y = posrot[:,1]
    s = npmovie.kalmanobj.speed[frames,0]
        
    cl = xy_trajectory(x,y,s,figure=figure, colormap=colormap)
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
    
    
    
    
def plot_movie_data(npmovie, show_wings=False, figure=None, legthresh=50, show_vision=True):
    calc_frame_of_landing(npmovie)
    frames = get_active_frames(npmovie)
    time = np.array(frames)*1/float(npmovie.fps)
    cl = xy_kalman(npmovie, figure=figure, frames=frames, colormap='gray')
    
    # get strobe image:
    strobe_img = strobe_from_npmovie(npmovie, interval=210, frames=frames)
    strobe_img = nim.rotate_image(strobe_img, np.array([[0,1],[-1,0]]))
    cl.ax0.imshow(strobe_img, pyplot.get_cmap('gray'), origin='lower')
    
    
    # axis parameters for subplots
    nxticks = 5
    nyticks = 3
    
    # subplot parameters
    n = 8
    h = (0.7-.05*(n-2))/float(n)
    
    subplots = []
    
    for s in range(n):
        ax = cl.fig.add_axes([0.47,0.1+(h+.05)*s,0.2,h])
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
    subplots[p].plot(time, npmovie.flycoord.speed[frames])
    subplots[p].set_xticks(xticks)
    default_yticks = subplots[p].get_yticks()
    yticks = np.linspace(default_yticks[0], default_yticks[-1], nyticks, endpoint=True).tolist()
    subplots[p].set_yticks(yticks)
    subplots[p].set_ylabel('speed')
    
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
        center = sa1_to_flydra_data_transform(npmovie.kalmanobj.positions[i], fix_sign=True)
        long_axis = sa1_to_flydra_img_transform(npmovie.kalmanobj.long_axis[i], fix_sign=False)*-1
                   
        factor =  npmovie.obj.axis_ratio[i][0]*4.
        factor = min(15., factor)
        
        dx = long_axis[1]*factor
        dy = long_axis[0]*factor
        
        legs = npmovie.kalmanobj.legs[i]
        
        #### get color from kmeans clusters ###
        speed = np.interp(npmovie.timestamps[i], npmovie.trajec.epoch_time, npmovie.trajec.speed)
        slipangle = npmovie.flycoord.slipangle[i]
        print npmovie.id, i
        dyaw = np.interp(npmovie.timestamps[i], npmovie.trajec.epoch_time[npmovie.sync2d3d.frames3d], npmovie.sync2d3d.smoothyaw[:,1])
        obs = np.nan_to_num(np.array([np.abs(dyaw), speed, np.abs(slipangle)]))
        cluster = kn.get_cluster_for_data(npmovie, obs)
        
        colormap = plt.get_cmap('jet')
        color = colormap((cluster) / float(npmovie.cluster_means.shape[0]))
        print 'color: ', color
        
        arrow = Arrow(center[0], center[1], dx, dy, width=1.0, color=color)
    
        cl.ax0.add_artist(arrow)
        
        if show_wings:
            # plot wing orientation vectors
            wingR = npmovie.kalmanobj.wingcenterR[i]
            if not np.isnan(wingR[0]):
                arrow = Arrow(center[0], center[1], wingR[1]-30, wingR[0]-30, width=1.0, color='b')
                cl.ax0.add_artist(arrow)
            
            wingL = npmovie.kalmanobj.wingcenterL[i]
            if not np.isnan(wingL[0]): 
                arrow = Arrow(center[0], center[1], wingL[1]-30, wingL[0]-30, width=1.0, color='b')
                cl.ax0.add_artist(arrow)
        
        
        i += interval
            
    #plt.show()
    try:
        if npmovie.vision_timeseries is not None:
            cl.vision_axes = [None for m in range(3)]
            for m in range(3):
                print 
                print m, len(npmovie.vision_timeseries)
                cl.vision_axes[m] = cl.fig.add_axes([0.75,0.1+0.25*m,0.2,.2])
                cl.vision_axes[m].imshow(npmovie.vision_timeseries[m], pyplot.get_cmap('gray'), origin='lower')
                cl.vision_axes[m].set_xticks([])
                time_raw = npmovie.timestamps[npmovie.sync2d3d.frames2d] - npmovie.timestamps[0]
                time = np.array([float(int(1000.*time_raw[i]))/1000. for i in range(len(time_raw))])

                cl.vision_axes[m].set_yticks(np.linspace(0,npmovie.vision_timeseries[m].shape[0], 5).tolist())
                cl.vision_axes[m].set_yticklabels(np.linspace(time[0], time[-1], 5).tolist())
                
                xlabel = 'horizontal slice of fly eye at ' + str( (m-1)*45 ) + ' deg'
                
                cl.vision_axes[m].set_xlabel(xlabel)
                cl.vision_axes[m].set_ylabel('time, sec')
    except:
        print 'no vision timeseries'    
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
            #center = npmovie.uframes[i].center
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



def pdf_movie_data(movie_dataset, scale = 10, movies = None):
    
    # Initialize:
    pp =  pdf.PdfPages('sa1_movies_2_test.pdf')
    f = 0
    
    if movies is None:
        movies = movie_dataset.keys()
    elif type(movies) is not list:
        movies = [movies]
       
        
    
    for key in movies:
        mnpmovie = movie_dataset[key]
        if mnpmovie.trajec is not None:
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
    return int(float(npmovie.fps)*timestamp)
        
#def align_sa1_flydra_fmin_func(dataset1, dataset2, index):
    
def interp_flydra_data_to_sa1(npmovie, data):
    flydra_raw_data = data
    flydra_raw_time = npmovie.trajec.fly_time
    flydra_interpolated_time = np.arange(flydra_raw_time[0], flydra_raw_time[-1], 1./float(npmovie.fps) )
    flydra_data = np.interp( flydra_interpolated_time, flydra_raw_time, flydra_raw_data )
    
    if return_time is False:
        return flydra_data
    else:
        return flydra_data, flydra_interpolated_time
    
def interp_sa1_data_to_flydra(npmovie, data, time, return_time=False):
    sa1_raw_data = data
    sa1_raw_time = time
    sa1_interpolated_time = np.arange(sa1_raw_time[0], sa1_raw_time[-1], 1./float(100.) )
    sa1_data = np.interp( sa1_interpolated_time, sa1_raw_time, sa1_raw_data )
    
    if return_time is False:
        return sa1_data
    else:
        return sa1_data, sa1_interpolated_time
    
def convolve(a1, a2):
    
    if len(a1) > len(a2):
        along = a1
        ashort = a2
    else:
        along = a2
        ashort = a1
    
    conv = np.zeros( [len(along)-len(ashort)] )
    for i in range(len(conv)):    
        conv[i] = np.sum(np.abs(along[i:i+len(ashort)] - ashort))
    
    return conv
    
    
def align_sa1_flydra(npmovie, plot=False):

    try:
        frames = get_active_frames(npmovie)
        first_frame_stamp = npmovie.timestamps[frames[0]]
        dt = npmovie.trajec.epoch_time - first_frame_stamp
        sa1_start_frame = np.argmin(np.abs(dt))
        npmovie.sa1_start_frame = sa1_start_frame
        return True

    except:
        print 'no timestamp data'
        return False
    
#    sa1_data = np.diff(interp_sa1_data_to_flydra(npmovie, npmovie.flycoord.dist_to_post[frames].T[0], npmovie.timestamps[frames]))
#    flydra_data = np.diff(npmovie.trajec.dist_to_stim_r + npmovie.trajec.stimulus.radius)
#    conv = convolve(sa1_data/sa1_data.max(), flydra_data/flydra_data.max())
    
    if plot:
        print flydra_data.shape, sa1_data.shape, conv.shape
        plt.plot(conv/np.max(conv))
        plt.plot(flydra_data/flydra_data.max())
        t = np.arange(np.argmax(conv), np.argmax(conv)+len(sa1_data), 1)
        plt.plot(t, sa1_data/sa1_data.max()) 
        plt.legend(['conv', 'flydra', 'sa1'])
    
#    index = np.argmax(conv)
#    if index < len(conv) and index > 0:
#        npmovie.sa1_start_index = np.argmax(conv)
#        return True
#    else:
#        npmovie.sa1_start_index = None
#        return False
        
        
def save_image_sequence(npmovie, frames, filename):
    # filename should be path + filename with initial number scheme + .extension desired, ie. jpg or png
    # the path should exist prior to calling the function
    imname = filename[filename.rfind('/')+1:]
    impath = filename[0:filename.rfind('/')+1]
    imname_no_extension = imname[0:imname.rfind('.')]
    imext = imname[imname.rfind('.'):]
    imbasename = imname_no_extension.rstrip('0123456789')
    imnumlen = len(imname_no_extension) - len(imbasename)

    i = -1
    for frame in frames:
        i += 1
        stri = str(i)
        while len(stri) < imnumlen:
            stri = '0' + stri        
        
        fname = impath + imbasename + stri + imext
        strobe = strobe_from_npmovie(npmovie, 1, frames=[frame, frame+1])
        strobe = nim.rotate_image(strobe, np.array([[0,1],[-1,0]]))
        plt.imsave(fname, strobe, origin='lower')
        
def get_frames_from_timestamps(npmovie, timestamps):
    frames = [timestamp2frame(npmovie, timestamp) for timestamp in timestamps]
    return frames    
    
def parse_path(filename):
    ''' splits a filename into the following:
    filename: '/home/floris/basename_001.jpg'
    path: '/home/floris/'
    basename: 'basename_'
    initnum: '001'
    ext: '.jpg'
    '''

    imname = filename[filename.rfind('/')+1:]
    impath = filename[0:filename.rfind('/')+1]
    imname_no_extension = imname[0:imname.rfind('.')]
    imext = imname[imname.rfind('.'):]
    imbasename = imname_no_extension.rstrip('0123456789')
    imnumlen = len(imname_no_extension) - len(imbasename)
    iminitnum = imname_no_extension[len(imbasename):]
    
    
    parsed = {}
    parsed.setdefault('path', impath)
    parsed.setdefault('basename', imbasename)
    parsed.setdefault('numlen', imnumlen)
    parsed.setdefault('initnum', iminitnum)
    parsed.setdefault('ext', imext)
    
    return parsed
    
def composite_movie():

    comp_im0 = '/home/floris/Desktop/composite/im_000.png'

    sa1_im0 = '/home/floris/Desktop/imseq/im_000.png'
    fsee_im0 = '/home/floris/src/fsee_utils/src/floris/fsee_sequence/ommatidia2image_reflines_00.png'
    
    sa1_fileparsing = parse_path(sa1_im0)
    fsee_fileparsing = parse_path(fsee_im0)
    comp_fileparsing = parse_path(comp_im0)
    
    pathparsings = [sa1_fileparsing, fsee_fileparsing]
    
    origins = [[0,0], [0,1024]]
    comp_imsize = [1024,1664]
    
    i = -1
    while 1:
        i += 1
        
        #if 1:
        try:
            images = []
            for pathparsing in pathparsings:
            
                initnum = pathparsing['initnum']
                n = i+int(initnum)
                strn = str(n)
                while len(strn) < pathparsing['numlen']:
                    strn = '0' + strn
                    
                filename = pathparsing['path'] + pathparsing['basename'] + strn + pathparsing['ext']
                im = plt.imread(filename)
                images.append(im)
        
            comp_rgb = nim.place_images(comp_imsize, images, origins)
            comp_filename = comp_fileparsing['path'] + comp_fileparsing['basename'] + strn + comp_fileparsing['ext']
            scipy.misc.imsave(comp_filename, comp_rgb)
        except:
            print 'cannot read files, iterations: ', i
            break
        
def flydra_to_sa1_transform(pt):
    m = np.array([[0,1],[-1,0]])
    if len(pt.shape) == 2:
        arr = np.zeros_like(pt)
        for i, p in enumerate(pt):
            arr[i] = np.dot(m,p)
        return arr
    else:
        return np.dot(m,pt)
    
def sa1_to_flydra_img_transform(pt, fix_sign=True):
    m = np.array([[0,-1],[1,0]])
    arr = sf_transform(pt, m, fix_sign=fix_sign)
    return arr
    
def sa1_to_flydra_data_transform(pt, fix_sign=False):
    m = np.array([[-1,0],[0,1]])
    arr = sf_transform(pt, m, fix_sign=fix_sign)
    return arr
    
def sf_transform(pt, m, fix_sign=False):
    if len(pt.shape) == 2:
        arr = np.zeros_like(pt)
        for i, p in enumerate(pt):
            arr[i] = np.dot(m,p)
            if fix_sign:    
                for j in range(2):
                    if arr[i,j] < 0:
                        arr[i,j] += 1024
        return arr
    else:
        arr = np.dot(m,pt)
        if fix_sign:    
                for j in range(2):
                    if arr[j] < 0:
                        arr[j] += 1024
        return arr
    
def normalize(arr):
    return (arr-arr.min()) / (arr-arr.min()).max()
    
def sync_2d_3d_data(npmovie, plot=False, res='lo'):
    
    if res == 'lo':
        is_sync_successful = align_sa1_flydra(npmovie, plot=False)
        if not is_sync_successful:
            return is_sync_successful
        else:
            npmovie.sync2d3d = pm.Object(npmovie)
            
            frames = get_active_frames(npmovie)
            interp_data, time = interp_sa1_data_to_flydra(npmovie, npmovie.flycoord.dist_to_post[frames].T[0], npmovie.timestamps[frames]-npmovie.timestamps[0], return_time=True)
            sync_frames_sa1 = get_frames_from_timestamps(npmovie, time)
            final_flydra_frame = np.min([npmovie.sa1_start_frame+len(sync_frames_sa1), len(npmovie.trajec.positions)])
            sync_frames_flydra = np.arange(npmovie.sa1_start_frame, final_flydra_frame)
            sync_frames_sa1 = sync_frames_sa1[0:len(sync_frames_flydra)]
            
            npmovie.sync2d3d.frames2d = sync_frames_sa1
            npmovie.sync2d3d.frames3d = sync_frames_flydra.tolist()
            
            npmovie.sync2d3d.pos2d = npmovie.kalmanobj.positions[npmovie.sync2d3d.frames2d,:]
            npmovie.sync2d3d.pos3d = npmovie.trajec.positions[npmovie.sync2d3d.frames3d,:]
            
            if plot is True:
                sa1_data = npmovie.flycoord.dist_to_post[npmovie.sync2d3d.frames2d].T[0]
                flydra_data = npmovie.trajec.dist_to_stim_r + npmovie.trajec.stimulus.radius
                plt.plot(npmovie.trajec.epoch_time, normalize(flydra_data))
                print npmovie.timestamps[npmovie.sync2d3d.frames2d].shape, sa1_data.shape, np.max(sa1_data)
                plt.plot(npmovie.timestamps[npmovie.sync2d3d.frames2d], normalize(sa1_data)) 
                plt.legend(['flydra', 'sa1'])
            return is_sync_successful

def calc_sa1_reconstructor(npmovie):
    pmat, residuals = camera_math.DLT(npmovie.sync2d3d.pos3d, npmovie.sync2d3d.pos2d, normalize = True)
    cal = flydra.reconstruct.SingleCameraCalibration(cam_id='sa1', Pmat=pmat, res=(1024,1024), scale_factor=1)
    npmovie.reconstructor = flydra.reconstruct.Reconstructor([cal])
    return npmovie.reconstructor
    
def set_sa1_reconstructor(npmovie, reconstructor):
    npmovie.reconstructor = reconstructor

    npmovie.sync2d3d = pm.Object(npmovie)
    frames = get_active_frames(npmovie)
    interp_data, time = interp_sa1_data_to_flydra(npmovie, npmovie.flycoord.dist_to_post[frames].T[0], npmovie.timestamps[frames], return_time=True)
    sync_frames_sa1 = get_frames_from_timestamps(npmovie, time)
    
    npmovie.sync2d3d.frames2d = sync_frames_sa1
    npmovie.sync2d3d.pos2d = npmovie.kalmanobj.positions[npmovie.sync2d3d.frames2d,:]

    try:    
        npmovie.sync2d3d.frames3d = sync_frames_flydra
        npmovie.sync2d3d.pos3d = npmovie.trajec.positions[npmovie.sync2d3d.frames3d,:]
    except:
        pass
    
    
def reproject_reconstruction(npmovie, estimated=False):
    plt.plot(npmovie.sync2d3d.pos2d[:,0], npmovie.sync2d3d.pos2d[:,1])
    reconstructed = np.zeros_like(npmovie.sync2d3d.pos2d)
    if estimated is False:
        data = npmovie.sync2d3d.pos3d
    else:
        data = npmovie.sync2d3d.pos3d_est
    for r in range(reconstructed.shape[0]):
        reconstructed[r,:] = npmovie.reconstructor.find2d('sa1', data[r])
    plt.plot(reconstructed[:,0], reconstructed[:,1], '*')
    return reconstructed
    
def calc_3d_estimate_from_2d(npmovie, z = 0):
    npmovie.sync2d3d.pos3d_est = np.zeros([npmovie.sync2d3d.pos2d.shape[0], 3])
    for r in range(npmovie.sync2d3d.pos3d_est.shape[0]):
        if type(z) is list:
            zval = z[r]
        else:
            zval = z
        npmovie.sync2d3d.pos3d_est[r] = npmovie.reconstructor.get_SingleCameraCalibration('sa1').get_3d_point_given_zval(npmovie.sync2d3d.pos2d[r], zval)
    
def calc_body_orientations(npmovie):
    def rotz(theta):
        ''' Returns a 3x3 rotation matrix corresponding to rotation around the *z* axis. '''
        return np.array([   [ np.cos(theta), -np.sin(theta), 0],
                            [ np.sin(theta), np.cos(theta), 0],
                            [0, 0, 1]])
                            
    npmovie.sync2d3d.long_axis = npmovie.kalmanobj.long_axis[npmovie.sync2d3d.frames2d,:]
    npmovie.sync2d3d.attitudes = [None for i in range(len(npmovie.sync2d3d.long_axis))]
    
    npmovie.sync2d3d.flydra_coordinates_longaxis = sa1_to_flydra_data_transform(npmovie.sync2d3d.long_axis)
    npmovie.sync2d3d.yaw = np.arctan2(npmovie.sync2d3d.flydra_coordinates_longaxis[:,1], npmovie.sync2d3d.flydra_coordinates_longaxis[:,0])
    for r in range(npmovie.sync2d3d.yaw.shape[0]):
        npmovie.sync2d3d.attitudes[r] = rotz(npmovie.sync2d3d.yaw[r])
    
    if 0: # this is an attempt at getting the right yaw out of the reprojection system... doesn't work
        npmovie.sync2d3d.long_axis = npmovie.kalmanobj.long_axis[npmovie.sync2d3d.frames2d,:]
        npmovie.sync2d3d.yaw = np.zeros([npmovie.sync2d3d.pos2d.shape[0]])
        npmovie.sync2d3d.attitudes = [None for i in range(len(npmovie.sync2d3d.yaw))]
        for r in range(npmovie.sync2d3d.yaw.shape[0]):
            center3d = npmovie.sync2d3d.pos3d_est[r]
            head2d = npmovie.sync2d3d.pos2d[r] + npmovie.sync2d3d.long_axis[r]*10
            head3d = npmovie.reconstructor.get_SingleCameraCalibration('sa1').get_example_3d_point_creating_image_point(head2d, 0.01)
            vector = head3d - center3d
            theta = np.arctan2(vector[1], vector[0])
            npmovie.sync2d3d.yaw[r] = theta
            npmovie.sync2d3d.attitudes[r] = rotz(theta)
        
def calc_function_for_dataset(movie_dataset, functions, args):
    # iterate through all npmovies for a list of functions
    for key,npmovie in movie_dataset.items():
        print 'processing: ', key
        for f, function in enumerate(functions):
            try:
                function(npmovie, args[f])
            except:
                function(npmovie)    
    
    
def prep_visual_timeseries_data(movie_dataset):
    for key,npmovie in movie_dataset.items():
        if npmovie.trajec is not None:
            print 'processing: ', key
            is_sync_successful = sync_2d_3d_data(npmovie)
            if is_sync_successful:
                calc_body_orientations(npmovie)
                try:
                    vision_timeseries = fsee_timeseries.calc_vision_timeseries(npmovie)
                except:
                    npmovie.vision_timeseries = None
            else:
                npmovie.vision_timeseries = None
        else:
            #print key
            npmovie.vision_timeseries = None
            
            
def get_sa1_timestamps_from_movie_dict(movie_dataset, movies):
    # movies should come from sa1_analysis.sa1_analysis()
    
    for k,npmovie in movie_dataset.items():
        npmovie.trigger_stamp = movies[k]['Trigger Stamp']
        npmovie.timestamps = movies[k]['Timestamps']
            
            
def smooth(data, Q, R, F=None, H=None, interpvals=0):
    os = data.shape[1] # observation size
    ss = os*2 # state size
    
    if F is None:
        if ss == 4:
            F = np.array([  [1,0,1,0], # process update
                            [0,1,0,1],
                            [0,0,1,0],
                            [0,0,0,1]],
                            dtype=np.float)
        elif ss == 2:
            F = np.array([   [1,1], # process update
                             [0,1]],
                             dtype=np.float)
                    
    if H is None:
        if ss == 4:
            H = np.array([   [1,0,0,0], # observation matrix
                             [0,1,0,0]],
                             dtype=np.float) 
                    
        elif ss == 2:         
            H = np.array([   [1,0] ], # observation matrix
                             dtype=np.float)
    
    data = np.nan_to_num(data)
    interpolated_data = np.zeros_like(data)
    
    for c in range(data.shape[1]):
        interpolated_data[:,c] = pm.interpolate(data[:,c], interpvals)
    y = interpolated_data
    initx = np.array([y[0,0],y[1,0]-y[0,0]],dtype=np.float)
    initV = 0*np.eye(ss)

    xsmooth,Vsmooth = adskalman.kalman_smoother(y,F,H,Q,R,initx,initV)

    return xsmooth, Vsmooth
            
            


