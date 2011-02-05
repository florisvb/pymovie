import sys

sys.path.append('/home/floris/src/adskalman')
sys.path.append('/home/floris/src/analysis')

import matplotlib.pyplot as plt
import numpy
from pyglet import media
import numpy as np
import pyglet
import numpyimgproc as nim
from scipy.ndimage.measurements import center_of_mass
from scipy.ndimage.morphology import binary_erosion
from scipy.ndimage.morphology import binary_dilation
import pickle as pickle
import copy

import numpy
import adskalman.adskalman as adskalman

THRESHHOLD = 15
THRESHRANGE = 15
WINGTHRESHRANGE = 15
ROI_RADIUS = 30
############################################################################################

def nextpow2(x):
    b = np.log2(x)
    p = np.ceil(b)
    return(p)

def calc_wingbeat_frequency(npmovie, framerange):
    
    arr = npmovie.kalmanobj.wingangleR[framerange[0]:framerange[1]]
    arr = arr.reshape(len(arr))
    L = len(arr)
    A = np.abs(np.fft.fft(arr, L))[1:int(L/2.)]    
    timestep = 1/npmovie.fps
    freq = np.fft.fftfreq(int(L), d=timestep)[1:int(L/2.)]    
    
    return freq[np.argmax(A)]

def load(filename):
    fname = (filename)
    fd = open( fname, mode='r')
    print 'loading data... from file'
    dataset = pickle.load(fd)
    return dataset

def save(movie, filename):
    fname = (filename)  
    fd = open( fname, mode='w' )
    pickle.dump(movie, fd)
    fd.close()
    return 1

def strobe_image(movie, interval=10, timerange=None, mode='darken'):
    # combine frames from a movie using lighten or darken modes to get a 'strobe-like' sequence
    # interval is the timestamp interval in movie-time (most likely 30 frames per second)

    if timerange is None:
        timerange = [0,movie.duration]
    
    strobe_image = None
    
    timestamp = timerange[0]
    while timestamp < timerange[1]:
        raw = movie.get_frame_at_timestamp(timestamp)
        if strobe_image is None:
            strobe_image = raw
        else:
            if mode == 'darken':
                strobe_image = nim.darken(strobe_image, raw)
            elif mode == 'lighten':
                strobe_image = nim.lighten(strobe_image, raw)
        timestamp += interval
        
    return strobe_image
    

###########################################################################################

#################################
#################################

def auto_adjust_uimg(npmovie):
    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe.uimg is not None:
            uframe.uimg_adj = nim.auto_adjust_levels(uframe.absdiff)
            
def calc_absdiff(npmovie):
    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe.uimg is not None:
            
            center = npmovie.obj.positions[i]
            x_lo = max( 0, int(center[0])-uframe.uimg.shape[0]/2 )
            x_hi = min( npmovie.background.shape[0], int(center[0])+uframe.uimg.shape[0]/2 )
            y_lo = max( 0, int(center[1])-uframe.uimg.shape[1]/2 )
            y_hi = min( npmovie.background.shape[0], int(center[1])+uframe.uimg.shape[1]/2 )
            bkgrd_uimg = npmovie.background[x_lo:x_hi, y_lo:y_hi]
            
            if uframe.uimg.shape != bkgrd_uimg.shape:
                uframe.absdiff = np.zeros_like(uframe.uimg)
            else:
                print 'good'
                uframe.absdiff = nim.absdiff(uframe.uimg, bkgrd_uimg)
            #except:
            #    
                
def smooth_axis_ratio(npmovie, radius=30):
    npmovie.kalmanobj.axis_ratio = np.zeros([len(npmovie.uframes), 1])
    for i in range(radius, len(npmovie.kalmanobj.axis_ratio)-radius):
        npmovie.kalmanobj.axis_ratio[i] = np.mean( npmovie.obj.axis_ratio[i-radius:i+radius+1] )
    npmovie.kalmanobj.axis_ratio /= np.max(npmovie.kalmanobj.axis_ratio)

def calc_polarpos(npmovie):
    
    center = npmovie.background.shape
    x0 = center[0]/2.
    y0 = center[1]/2.
    
    npmovie.kalmanobj.polarpos = np.zeros_like(npmovie.kalmanobj.positions)
    
    for i in range(len(npmovie.kalmanobj.polarpos)):
        x = npmovie.kalmanobj.positions[i][0]-x0
        y = npmovie.kalmanobj.positions[i][1]-y0
        r = np.sqrt((x)**2 + (y)**2)
        
        vec_to_post = np.array([x0-x, y0-y]) / np.linalg.norm(np.array([x0-x, y0-y]))
        body_angle = npmovie.kalmanobj.long_axis[i]
        
        theta = np.arccos(np.cross(vec_to_post, body_angle))
        npmovie.kalmanobj.polarpos[i] = [r, theta]
            
def find_object(img, thresholds, sizerange, dist_thresh, erode=False, check_centers=True):
    
            body = nim.threshold(img, thresholds[0], thresholds[1])
            
            if erode is not False:
                for i in range(erode):
                    body = binary_erosion(body)
                    
            if check_centers is False:
                blobs = nim.find_blobs(body, sizerange=sizerange, aslist=False)
            else:
                blobs = nim.find_blobs(body, sizerange=sizerange, aslist=True)
            body = blobs
            
            if check_centers:
                centers = nim.center_of_blob(blobs)
                dist = []
                for center in centers:
                    diff = np.linalg.norm( center - np.array(img.shape)/2. )
                    dist.append(diff)
                body = np.zeros_like(img)
                for j, d in enumerate(dist):
                    if d < dist_thresh:
                        body += blobs[j]
                                                
            if body.max() > 1:
                body /= body.max()
            
            if body is None:
                body = np.zeros_like(img)
            
            body = np.array(body*255, dtype=np.uint8)
            
            
            return body
            
                        
def segment_fly(npmovie):

    npmovie.obj.legs = np.zeros(len(npmovie.uframes))

    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe.uimg is not None:
            body = find_object(uframe.uimg_adj, [150,253], [100,400], 10, check_centers=False)
            uframe.body = copy.copy(body)
            
            raw_wing_img = copy.copy(uframe.uimg_adj)
            body = binary_dilation(body)
            raw_wing_img[body>0] = 255
            wings = find_object(raw_wing_img, [30,140], [15,400], 17, erode=1, check_centers=False)
            uframe.wings = copy.copy(wings)
            
            if 1:
                raw_leg_img = copy.copy(uframe.uimg_adj)
                wings = binary_dilation(wings)
                wings = binary_dilation(wings)
                wings = binary_dilation(wings)
                body = binary_dilation(body)
                body = binary_dilation(body)
                raw_leg_img[body>0] = 255
                raw_leg_img[wings>0] = 255
                uframe.raw_leg_img = raw_leg_img
                
                legs = find_object(raw_leg_img, [30,120], [6,100], 15, erode=False, check_centers=True)
                uframe.legs = copy.copy(legs)
                
                npmovie.obj.legs[i] = uframe.legs.sum()
                
                
def smooth_legs(npmovie, clipping_radius=7, smoothing_radius=200):
    npmovie.kalmanobj.legs = np.zeros([len(npmovie.uframes), 1])
    for i in range(clipping_radius, len(npmovie.kalmanobj.legs)-clipping_radius):
        npmovie.kalmanobj.legs[i] = np.min( npmovie.obj.legs[i-clipping_radius:i+clipping_radius+1] )
    for i in range(smoothing_radius, len(npmovie.kalmanobj.legs)-smoothing_radius):
        npmovie.kalmanobj.legs[i] = np.mean( npmovie.kalmanobj.legs[i-smoothing_radius:i+smoothing_radius+1] )
    #npmovie.kalmanobj.legs /= np.max(npmovie.kalmanobj.legs)



#################################
#################################




############################################################################################
class Object:
    def __init__(self, npmovie):
        self.positions = np.zeros([len(npmovie.uframes), 2])
        self.velocities = np.zeros([len(npmovie.uframes), 2])
        self.long_axis = np.zeros([len(npmovie.uframes), 2])
        self.axis_ratio = np.zeros([len(npmovie.uframes), 1])
        self.speed = np.zeros([len(npmovie.uframes), 1])

def calc_obj_motion(npmovie):
    npmovie.kalmanobj = Object(npmovie)

    ss = 4 # state size
    os = 2 # observation size
    F = numpy.array([[1,0,1,0], # process update
                     [0,1,0,1],
                     [0,0,1,0],
                     [0,0,0,1]],
                    dtype=numpy.float)
    H = numpy.array([[1,0,0,0], # observation matrix
                     [0,1,0,0]],
                    dtype=numpy.float)
    Q = 0.0001*numpy.eye(ss) # process noise
    R = 100*numpy.eye(os) # observation noise
    
    raw_x = np.nan_to_num(npmovie.obj.positions[:,0])
    indices = np.nonzero(raw_x)[0].tolist()
    raw_x = raw_x[indices]
    
    raw_y = np.nan_to_num(npmovie.obj.positions[:,1])
    raw_y = raw_y[indices]
    
    y = np.vstack( (raw_x, raw_y) ).T
    print y.shape, raw_x.shape, raw_y.shape
    initx = numpy.array([y[0,0], y[0,1],0,0],dtype=numpy.float)
    initV = 0*numpy.eye(ss)

    xsmooth,Vsmooth = adskalman.kalman_smoother(y,F,H,Q,R,initx,initV)

    npmovie.kalmanobj.positions[indices] = xsmooth[:,0:2]
    npmovie.kalmanobj.velocities[indices] = xsmooth[:,2:4]
    npmovie.kalmanobj.timestamps = [npmovie.uframes[i].timestamp for i in indices]
    npmovie.kalmanobj.indices = indices
    npmovie.kalmanobj.errors = Vsmooth
    
    npmovie.kalmanobj.speed = np.zeros([len(npmovie.uframes), 1])
    for i, v in enumerate(npmovie.kalmanobj.velocities):
        npmovie.kalmanobj.speed[i] = np.linalg.norm(v)
    
    # need to fix/smooth missing angles
    for i in range(len(npmovie.kalmanobj.indices)):
        frame = indices[i]
        npmovie.kalmanobj.long_axis[frame] = npmovie.obj.long_axis[frame] / np.linalg.norm(npmovie.obj.long_axis[frame])
    for i in range(1,len(npmovie.kalmanobj.indices)):
        frame = indices[i]
        if npmovie.kalmanobj.long_axis[frame][0] == 1:
            
            future_axis = 1
            future_frame = frame
            while future_axis == 1:
                future_frame += 1
                if future_frame > len(npmovie.kalmanobj.indices):
                    future_frame = frame
                    break
                future_axis = npmovie.kalmanobj.long_axis[future_frame][0]
            
            delta_axis = (npmovie.kalmanobj.long_axis[future_frame] - npmovie.kalmanobj.long_axis[frame-1]) / float(future_frame- (frame-1) )
            
            npmovie.kalmanobj.long_axis[frame] = npmovie.kalmanobj.long_axis[frame-1] + delta_axis
            
        
    
    # fix angle orientation:
    if 1:
        switching_threshold = 0.004
        for i in range(len(npmovie.kalmanobj.indices)):
            
            # find the value of orienting the body with velocity
            value_new = np.dot(npmovie.kalmanobj.velocities[frame], npmovie.obj.long_axis[frame])
            if value_new < 0:
                npmovie.obj.long_axis[frame] *= -1
                value_new = np.dot(npmovie.kalmanobj.velocities[frame], npmovie.obj.long_axis[frame])
            
            value_prev = np.dot(npmovie.kalmanobj.velocities[frame], npmovie.kalmanobj.long_axis[frame-1])

            print 'values: ', value_new, value_prev, switching_threshold

            if value_new > value_prev-switching_threshold:
                npmovie.kalmanobj.long_axis[frame] = npmovie.obj.long_axis[frame]
            else:
                npmovie.kalmanobj.long_axis[frame] = npmovie.kalmanobj.long_axis[frame-1]
            
    return npmovie

        
def calc_obj_pos(npmovie):

    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe.uimg is not None:
            npmovie.obj.positions[i,:] = uframe.center
            npmovie.obj.long_axis[i,:], npmovie.obj.axis_ratio[i] = nim.fit_ellipse_cov(uframe.body)
                
    return npmovie
################################################################################################################

class Movie:
    def __init__(self, filename=None):
        self.mask = None
        self.filename = filename
        if filename is not None:
            self.load_source(filename)
    
    def load_source(self, filename):
        print 'loading source'
        try:
            self.source = media.load(filename)
            self.width = int(self.source.video_format.width)
            self.height = int(self.source.video_format.height)
            self.duration = self.source.duration
            print 'loaded movie.. height: ', self.height, ' | width: ', self.width
        except:
            print 'failed to open movie!'
            ValueError('failed to open movie!')
            
    def get_next_frame(self):
        imdata = self.source.get_next_video_frame()
        a = numpy.frombuffer(imdata.data, numpy.uint8)
        a.shape = (imdata.height, imdata.width, 3)
        del(imdata)
        raw = (a[:,:,0]+a[:,:,1]+a[:,:,2])/3 #convert to rawchrome
        raw = copy.copy(raw)
        del(a)
        if self.mask is not None:
            raw += np.ones_like(raw)
            raw *= self.mask
            # delete rows that are all zeros:
            nz = (raw == 0).sum(1)
            raw = raw[nz == 0, :]
            # delete columns that are all zeros:
            nz = (raw == 0).sum(0)
            raw = raw[:, nz == 0]
            raw -= np.ones_like(raw)
        return raw
        
    def get_frame_at_timestamp(self, timestamp):
        self.seek(timestamp)
        raw = self.get_next_frame()
        return raw
            
    def seek(self, timestamp=0):
        self.source._seek(timestamp)

class npMovie:
    def __init__(self):
        self.frames = []
        self.uframes = []
        
        self.width = None
        self.height = None
        
        self.tracking_mask = None
        
        self.dynamic_tracking_mask = False
        self.mask_center = [0,0]
        self.mask_radius = 50
        
        self.blob_size_range = [50,400]
        
        self.fps = None
        
    def calc_background(self, movie):
        
        first_frame = movie.get_frame_at_timestamp(0)
        last_frame = movie.get_frame_at_timestamp(movie.duration-1)
        self.background = nim.lighten(first_frame, last_frame)
        self.height, self.width = self.background.shape
        
    def background_subtraction(self, raw, save_raw=False):
    
        mask_center_0 = int(self.mask_center[0])
        mask_center_1 = int(self.mask_center[1])
        mask_0_lo = max(0, mask_center_0 - self.mask_radius)
        mask_0_hi = min(self.width, mask_center_0 + self.mask_radius)
        mask_1_lo = max(0, mask_center_1 - self.mask_radius)
        mask_1_hi = min(self.width, mask_center_1 + self.mask_radius)
        if self.tracking_mask is not None:
            tracking_mask = copy.copy(self.tracking_mask)
            tracking_mask[mask_0_lo:mask_0_hi, mask_1_lo:mask_1_hi] = 1
        
        if self.dynamic_tracking_mask is True:
            # TODO: currently dynamic tracking and such only works for square format cameras 
            print 'dynamic tracking'
            masked_img = raw[mask_0_lo:mask_0_hi, mask_1_lo:mask_1_hi]
            masked_background = self.background[mask_0_lo:mask_0_hi, mask_1_lo:mask_1_hi]
        else:
            masked_img = raw
            masked_background = self.background

        if masked_img.shape[0] < 100 or masked_img.shape[1] < 100:
            print 'no uframe'
            self.dynamic_tracking_mask = False
            uframe = uFrame()
            return uframe, None        
            
        absdiff = nim.absdiff(masked_img, masked_background)
        if self.dynamic_tracking_mask is False and self.tracking_mask is not None:
            absdiff *= tracking_mask
            
        absdiff = nim.auto_adjust_levels(absdiff)
        #print 'shape: ', np.shape(absdiff)
        #threshold = max( 10, absdiff.max() - THRESHRANGE )
        #print 'dynamic threshold: ', threshold 
        
        #diffthresh = nim.threshold(absdiff, threshold)
        diffthresh = nim.threshold(absdiff, 50, threshold_hi=255)
        blobs = nim.find_blobs(diffthresh, self.blob_size_range)
        
        if blobs is None:
            self.dynamic_tracking_mask = False
            masked_img = raw
            masked_background = self.background
            
            absdiff = nim.absdiff(masked_img, masked_background)
            diffthresh = nim.threshold(absdiff, threshold)
            if self.dynamic_tracking_mask is False and self.tracking_mask is not None:
                diffthresh *= tracking_mask
            blobs = nim.find_blobs(diffthresh, self.blob_size_range)
        
        
        if blobs is None:
            uframe = uFrame()
            if save_raw is False:
                print 'no uframe'
                self.dynamic_tracking_mask = False
                return uframe, None
            else:
                frame = Frame(raw, absdiff, diffthresh)
                return uframe, frame
            
        nblobs = blobs.max()
        if nblobs > 1:
            blobs = nim.find_biggest_blob(blobs)
            
        center = nim.center_of_blob(blobs)
        
        if np.isnan(center)[0] or np.isnan(center)[1]:
            uframe = uFrame()
            print 'no uframe, NaN center!'
            self.dynamic_tracking_mask = False
            return uframe, None
        if center[0] < 1 or center[1] < 1:
            uframe = uFrame()
            print 'no uframe, NaN center!'
            self.dynamic_tracking_mask = False
            return uframe, None
        
        print 'center found'
        if 1:
            limlo_x = max( int(center[0])-ROI_RADIUS, 0 )
            limlo_y = max( int(center[1])-ROI_RADIUS, 0 )
            
            limhi_x = min( int(center[0])+ROI_RADIUS, masked_img.shape[0] )
            limhi_y = min( int(center[1])+ROI_RADIUS, masked_img.shape[1] )
            
            # TODO: right now object entering or leaving frame doesnt work perfectly
            uimg = masked_img[limlo_x:limhi_x, limlo_y:limhi_y]
            uimg = copy.copy(uimg)
            uimg_absdiff = absdiff[limlo_x:limhi_x, limlo_y:limhi_y]
            uimg_absdiff = copy.copy(uimg_absdiff)
            uimg_diffthresh = blobs[limlo_x:limhi_x, limlo_y:limhi_y]
            uimg_diffthresh = copy.copy(uimg_diffthresh)
            
            if self.dynamic_tracking_mask is True:
                tmp = np.array([masked_img.shape[0]/2., masked_img.shape[1]/2.])
                center = np.array(self.mask_center) + np.array(center) - tmp
                uimg_indices = [mask_0_lo+limlo_x, mask_0_lo+limhi_x, mask_1_lo+limlo_y, mask_1_lo+limhi_y]
            else:
                uimg_indices = [limlo_x, limhi_x, limlo_y, limhi_y]
            self.dynamic_tracking_mask = True
            self.mask_center = center
            
            uframe = uFrame(center, uimg, uimg_absdiff, uimg_diffthresh) 
            uframe.blobsize = blobs.sum()
            
            uframe.indices = uimg_indices
                
            if save_raw is False:
                del(raw)
                del(absdiff)
                del(diffthresh)
                del(blobs)
                return uframe, None
                
            frame = Frame(masked_img, absdiff, diffthresh)
        return uframe, frame    
        
    def load_frames(self, movie, timerange=None, save_raw=False):
        print 'loading frames.. '
        
        if timerange is None:
            timerange = [0,movie.duration]
        if timerange[0] == 0:
            timerange[0] = 0.1 # hack
        movie.seek(timerange[0])
        
        print 'timerange: ', timerange[0], ' to ', timerange[1]
            
        fr = -1
        while True:
        
            timestamp = movie.source.get_next_video_timestamp()
            if timestamp > timerange[1]:
                break
            try:
                raw = movie.get_next_frame()    
            except:
                break
            if 1:       
                uframe,frame = self.background_subtraction(raw, save_raw=save_raw)
                uframe.timestamp = timestamp
                del(raw)
                
                if save_raw is True:
                    frame.timestamp = timestamp
                    self.frames.append(frame)
                else:
                    del(frame)
                
                self.uframes.append(uframe)
            fr += 1
            print 'loading...', movie.filename, ' frame num: ', fr, '  ::  ', 'timestamp: ', timestamp, ' :: ', timestamp/movie.duration*100, '%'
            
class Frame:
    def __init__(self, raw, absdiff=None, diffthresh=None):
        self.raw = raw
        self.size = np.shape(self.raw)
        self.width = self.size[1]
        self.height = self.size[0]
        self.absdiff = absdiff
        self.diffthresh = diffthresh
    def show(self, var='raw'):
        if var == 'raw':
            plt.imshow(self.raw)
        elif var == 'absdiff':
            plt.imshow(self.absdiff)
        elif var == 'diffthresh':
            plt.imshow(self.diffthresh)  
        else:
            print 'please use a valid var'
            return
        plt.show()
        return
        
class uFrame:
    def __init__(self, center=None, uimg=None, absdiff=None, diffthresh=None):
        self.uimg = uimg
        self.center = center
        self.absdiff = absdiff
        self.diffthresh = diffthresh
        self.wingimg = None
        self.timestamp = None
        self.flydraframe = None
    def show(self, var='uimg'):
        if var == 'uimg':
            plt.imshow(self.uimg)
        elif var == 'absdiff':
            plt.imshow(self.absdiff)
        elif var == 'diffthresh':
            plt.imshow(self.diffthresh)
        elif var == 'wingimg':
            plt.imshow(self.wingimg)  
        elif var == 'wingimg2':
            plt.imshow(self.wingimg2)
        elif var == 'flysegs':
            plt.imshow(self.flysegs)
        else:
            print 'please use a valid var'
            return
        plt.show()
        return        
        
class MiniNPM:
    def __init__(self, npmovie):
        self.background = copy.copy(npmovie.background)
        self.obj = copy.copy(npmovie.obj)
        self.kalmanobj = copy.copy(npmovie.kalmanobj)
        self.fps = copy.copy(npmovie.fps)
        self.trajec = copy.copy(npmovie.trajec)
        try:
            self.dataset_id = copy.copy(npmovie.dataset_id)
        except:
            self.dataset_id = None
        self.epochtime = copy.copy(npmovie.epochtime)
        
        self.uframes = [uFrame() for i in range(len(npmovie.uframes))]
        
        for i, uframe in enumerate(self.uframes):
            uframe.uimg = copy.copy(npmovie.uframes[i].uimg) 
            uframe.center = copy.copy(npmovie.uframes[i].center)
            try:
                self.uframes[i].flydraframe = copy.copy(npmovie.uframes[i].flydraframe)
                uframe.indices = copy.copy(npmovie.uframes[i].indices)
                #print i, uframe.flydraframe
            except:
                self.uframes[i].flydraframe = self.uframes[i-1].flydraframe
                uframe.indices = None
    
                
if __name__ == '__main__':

    filename =   '/media/SA1_movies_2/sa1_movies/20101028_C001H001S0008.avi'
    
    movie = Movie(filename)
    npmovie = npMovie()
    
    # make a mask
    mask = np.ones([movie.height, movie.width])
    mask[0:movie.height-movie.width, :] = 0
    movie.mask = mask
    
    # make a tracking mask -- 0 where no tracking desired
    tracking_mask = np.abs(nim.plot_circle(1024,1024,[512,512], 150)-1)
    npmovie.tracking_mask = tracking_mask
    
    npmovie.calc_background(movie)
    npmovie.load_frames(movie)
    
    npmovie.obj = Object(npmovie)
    calc_obj_pos(npmovie)
    calc_obj_motion(npmovie, 5000.)
    save(npmovie, 'sa1_movie_obj_flyby_20101028_C001H001S0008')
    
    
