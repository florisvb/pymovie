import sys

sys.path.append('/home/floris/src/adskalman')

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
WINGTHRESHRANGE = 17
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

def find_wings(npmovie):
    
    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe.uimg is not None:
            if uframe.uimg.shape[0] == uframe.uimg.shape[1]:
                try:
                    uimg = uframe.uimg
                    body = binary_dilation(uframe.diffthresh)
                    body = binary_dilation(body)
                    uimg_no_body = uframe.uimg + body*255
                    
                    threshold = uimg_no_body.min()+WINGTHRESHRANGE
                    
                    uimg_no_body_thresh = uimg_no_body<threshold
                    wingimg = binary_erosion(uimg_no_body_thresh)
                    uframe.wingimg = wingimg
                except:
                    pass         
            else:
                uframe.wingimg = None 
                
def split_wings_bodyaxis(npmovie):
    
    for i in range(len(npmovie.uframes)):
        #i = 11
        uframe = npmovie.uframes[i]
        if uframe.wingimg is not None:
            wingimg = copy.copy(np.array(npmovie.uframes[i].wingimg, dtype=int))  
            #wingimg = np.ones_like(wingimg)
            #wingimg = copy.copy(np.array(npmovie.uframes[i].diffthresh, dtype=int))  
            long_axis = npmovie.kalmanobj.long_axis[i]
            
            
            x1 = np.array([30, 30])
            x2 = np.array([30+long_axis[1], 30+long_axis[0]])
            
            m = long_axis[1]/long_axis[0]
            b = x2[1]-m*x2[0]
            
            def above_body_axis(x,y):
                pt_on_line = m*x+b
                if y > pt_on_line:
                    return True
                else:
                    return False
                    
            for x in range(wingimg.shape[0]):
                for y in range(wingimg.shape[1]):
                    is_above_axis = above_body_axis(x,y)
                    if is_above_axis:
                        wingimg[x,y] *= -1
            
            npmovie.uframes[i].wingimg2 = np.array((wingimg*127)+127, dtype=np.uint8)
            npmovie.uframes[i].flysegs = np.array((wingimg*127)+127 + np.array(npmovie.uframes[i].diffthresh,dtype=np.uint8)*50, dtype=np.uint8)
            npmovie.uframes[i].wingimgR = np.array( (wingimg>0) , dtype=np.uint8) *255
            npmovie.uframes[i].wingimgL = np.array( (wingimg<0) , dtype=np.uint8) *255
            
def split_wings_blobs(npmovie):
    
    for i in range(len(npmovie.uframes)):
        #i = 11
        uframe = npmovie.uframes[i]
        if uframe.wingimg is not None:
            wingimg = copy.copy(np.array(npmovie.uframes[i].wingimg, dtype=int))  
            
            #wingimg = np.ones_like(wingimg)
            #wingimg = copy.copy(np.array(npmovie.uframes[i].diffthresh, dtype=int))  
            long_axis = npmovie.kalmanobj.long_axis[i]
            x1 = np.array([30, 30])
            x2 = np.array([30+long_axis[1], 30+long_axis[0]])
            m = long_axis[1]/long_axis[0]
            b = x2[1]-m*x2[0]
            
            def above_body_axis(x,y):
                pt_on_line = m*x+b
                if y > pt_on_line:
                    return True
                else:
                    return False
                    
            def get_mean_axis(w1):
                for x in range(w1.shape[0]):
                    for y in range(w1.shape[1]):
                        is_above_axis = above_body_axis(x,y)
                        if is_above_axis:
                            w1[x,y] *= -1
                w1_avg = np.mean(w1)
                return w1_avg
                
            # find biggest blob, check to see if its center is near where the center of a wing was last round. if so, it's the new position of that wing. 
                
            w1, w2 = nim.find_2_biggest_blobs(wingimg)
            
            if w1 is not None:
                c1 = nim.center_of_blob(w1)
                c1_above_axis = above_body_axis(c1[0],c1[1])
                w1_avg = get_mean_axis(w1)
            else:
                w1_avg = 0
                w1 = np.zeros_like(wingimg)
                
            if w2 is not None:
                c2 = nim.center_of_blob(w2)
                c2_above_axis = above_body_axis(c2[0],c2[1])
                w2_avg = get_mean_axis(w2)
            else:
                w2_avg = 0
                w2 = np.zeros_like(wingimg)
            
            if c1_above_axis:
                npmovie.uframes[i].wingimgR = np.array( w1 , dtype=np.uint8) 
                npmovie.uframes[i].wingimgL = np.array( w2 , dtype=np.uint8) 
            elif c2_above_axis:
                npmovie.uframes[i].wingimgR = np.array( w2 , dtype=np.uint8) 
                npmovie.uframes[i].wingimgL = np.array( w1 , dtype=np.uint8) 
            else:
                npmovie.uframes[i].wingimgR = np.array( np.zeros_like(wingimg) , dtype=np.uint8) 
                npmovie.uframes[i].wingimgL = np.array( np.zeros_like(wingimg) , dtype=np.uint8) 
                
            npmovie.uframes[i].wingimg2 = np.array((wingimg*127)+127, dtype=np.uint8)
            npmovie.uframes[i].flysegs = np.array(  npmovie.uframes[i].wingimgR*127 +
                                                    npmovie.uframes[i].wingimgL*255 +
                                                    np.array(npmovie.uframes[i].diffthresh,dtype=np.uint8)*50, dtype=np.uint8)
            

def split_wings_blobs(npmovie):
    
    for i in range(len(npmovie.uframes)):
        #i = 11
        uframe = npmovie.uframes[i]
        if uframe.wingimg is not None:
            wingimg = copy.copy(np.array(npmovie.uframes[i].wingimg, dtype=int))  
            
            #wingimg = np.ones_like(wingimg)
            #wingimg = copy.copy(np.array(npmovie.uframes[i].diffthresh, dtype=int))  
            long_axis = npmovie.kalmanobj.long_axis[i]
            x1 = np.array([30, 30])
            x2 = np.array([30+long_axis[1], 30+long_axis[0]])
            m = long_axis[1]/long_axis[0]
            b = x2[1]-m*x2[0]
            
            def above_body_axis(x,y):
                pt_on_line = m*x+b
                if y > pt_on_line:
                    return True
                else:
                    return False
                    
            def get_mean_axis(w1):
                for x in range(w1.shape[0]):
                    for y in range(w1.shape[1]):
                        is_above_axis = above_body_axis(x,y)
                        if is_above_axis:
                            w1[x,y] *= -1
                w1_avg = np.mean(w1)
                return w1_avg
                
            # find biggest blob, check to see if its center is near where the center of a wing was last round. if so, it's the new position of that wing. 
                
            w1, w2 = nim.find_2_biggest_blobs(wingimg)
            
            if w1 is not None:
                c1 = nim.center_of_blob(w1)
                c1_above_axis = above_body_axis(c1[0],c1[1])
                w1_avg = get_mean_axis(w1)
            else:
                w1_avg = 0
                w1 = np.zeros_like(wingimg)
                
            if w2 is not None:
                c2 = nim.center_of_blob(w2)
                c2_above_axis = above_body_axis(c2[0],c2[1])
                w2_avg = get_mean_axis(w2)
            else:
                w2_avg = 0
                w2 = np.zeros_like(wingimg)
            
            if c1_above_axis:
                npmovie.uframes[i].wingimgR = np.array( w1 , dtype=np.uint8) 
                npmovie.uframes[i].wingimgL = np.array( w2 , dtype=np.uint8) 
            elif c2_above_axis:
                npmovie.uframes[i].wingimgR = np.array( w2 , dtype=np.uint8) 
                npmovie.uframes[i].wingimgL = np.array( w1 , dtype=np.uint8) 
            else:
                npmovie.uframes[i].wingimgR = np.array( np.zeros_like(wingimg) , dtype=np.uint8) 
                npmovie.uframes[i].wingimgL = np.array( np.zeros_like(wingimg) , dtype=np.uint8) 
                
            npmovie.uframes[i].wingimg2 = np.array((wingimg*127)+127, dtype=np.uint8)
            npmovie.uframes[i].flysegs = np.array(  npmovie.uframes[i].wingimgR*127 +
                                                    npmovie.uframes[i].wingimgL*255 +
                                                    np.array(npmovie.uframes[i].diffthresh,dtype=np.uint8)*50, dtype=np.uint8)
                                                    
                                                        
            
        
def calc_wing_motion(npmovie):
    npmovie.obj.wingcenterR = np.zeros([len(npmovie.uframes), 2])
    #npmovie.obj.wingaxisR = np.zeros([len(npmovie.uframes), 2])

    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe.wingimg is not None:
            center, long_axis, ratio = nim.fit_ellipse_cov(uframe.wingimgR, recenter=True)
            npmovie.obj.wingcenterR[i] = center
            
            
def add_wings(npmovie, mode='lighten'):
    
    wingsum = None
    
    for uframe in npmovie.uframes:
        if uframe is not None:
            if wingsum is None:
                wingsum = uframe.wingimg
            else:
                if mode == 'darken':
                    wingsum = nim.darken(wingsum, uframe.wingimg)
                elif mode == 'lighten':
                    wingsum = nim.lighten(wingsum, uframe.wingimg)
                    
    return wingsum
    
def calc_wing_angle(npmovie):
    npmovie.obj.wingangleR = np.zeros([len(npmovie.uframes), 1])
    
    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe.wingimg is not None:
            wingvector = npmovie.kalmanobj.wingaxisR[i] 
            bodyvector = npmovie.kalmanobj.long_axis[i] 
            dotprod = np.dot(wingvector, bodyvector)
            cosangle = dotprod / np.linalg.norm(wingvector)*np.linalg.norm(bodyvector)
            angle = np.arccos(cosangle)
            if np.isnan(angle):
                try:
                    angle = npmovie.obj.wingangleR[i-1]
                except:
                    angle = 0.
            npmovie.obj.wingangleR[i] = angle
            
def interp_nan_values(npmovie):
    warr = npmovie.obj.wingcenterR
    npmovie.kalmanobj.wingcenterR = np.zeros([len(npmovie.uframes), 2])
    npmovie.kalmanobj.wingaxisR = np.zeros([len(npmovie.uframes), 2])
    
    for i, w in enumerate(warr):
        
        if np.isnan(w[0]):
            prev = npmovie.kalmanobj.wingcenterR[i-1]
            next = np.array([np.nan, np.nan])
            n = 1
            while np.isnan(next[0]):
                next = warr[i+n]
                n += 1
            diff = next-prev
            inc = diff/float(n)
            new = prev+inc
            npmovie.kalmanobj.wingcenterR[i] = new
        else:
            npmovie.kalmanobj.wingcenterR[i] = warr[i]
            
    vec = npmovie.kalmanobj.wingcenterR - npmovie.uframes[i].uimg.shape[0]
    for i, v in enumerate(vec):
        npmovie.kalmanobj.wingaxisR[i] = np.nan_to_num( vec[i] / np.linalg.norm(vec[i]) )
        
                
def kalmanize_wing_angle(npmovie):

    npmovie.kalmanobj.wingangleR = np.zeros([len(npmovie.uframes), 1])

    ss = 2 # state size
    os = 1 # observation size
    F = numpy.array([[1,0],
                     [0,1]],
                    dtype=numpy.float)
    H = numpy.array([[1,0]],
                    dtype=numpy.float)
    Q = 0.0001*numpy.eye(ss) # process noise
    R = 0.001*numpy.eye(os) # observation noise
    
    raw_x = npmovie.obj.wingangleR
    indices = np.nonzero(raw_x)[0].tolist()
    raw_x = raw_x[indices]
    
    initx = numpy.array([raw_x[0][0],0],dtype=numpy.float)
    initV = 0*numpy.eye(ss)

    xsmooth,Vsmooth = adskalman.kalman_smoother(raw_x,F,H,Q,R,initx,initV)
    
    npmovie.kalmanobj.wingangleR[indices] = xsmooth[:,0].reshape(xsmooth.shape[0], 1)

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
    
    initx = numpy.array([y[0,0], y[0,1],0,0],dtype=numpy.float)
    initV = 0*numpy.eye(ss)

    xsmooth,Vsmooth = adskalman.kalman_smoother(y,F,H,Q,R,initx,initV)

    npmovie.kalmanobj.positions[indices] = xsmooth[:,0:2]
    npmovie.kalmanobj.velocities[indices] = xsmooth[:,2:4]
    npmovie.kalmanobj.timestamps = [npmovie.uframes[i].timestamp for i in indices]
    npmovie.kalmanobj.indices = indices
    npmovie.kalmanobj.speed = np.linalg.norm(npmovie.kalmanobj.velocities[i,:]) 
    npmovie.kalmanobj.errors = Vsmooth
    
    # fix angles:
    dot_arr = np.zeros([len(npmovie.uframes), 1])
    
    if 1:
        for i in range(len(npmovie.kalmanobj.indices)):
            frame = indices[i]
            npmovie.kalmanobj.long_axis[frame] = npmovie.obj.long_axis[frame] / np.linalg.norm(npmovie.obj.long_axis[frame])
            dot = np.dot(npmovie.kalmanobj.velocities[frame], npmovie.obj.long_axis[frame])
            dot_arr[frame] = np.sign(dot)
        mode_dot = np.sign(np.mean(dot_arr))
        if mode_dot < 0:
            npmovie.kalmanobj.long_axis *= -1
            
    return npmovie

        
def calc_obj_pos(npmovie):

    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe is not None:
            
            npmovie.obj.positions[i,:] = uframe.center
            npmovie.obj.long_axis[i,:], npmovie.obj.axis_ratio[i] = nim.fit_ellipse_cov(uframe.diffthresh)
                
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
        except:
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
            tracking_mask = self.tracking_mask
            tracking_mask[mask_0_lo:mask_0_hi, mask_1_lo:mask_1_hi] = 1
        
        if self.dynamic_tracking_mask is True:
            # TODO: currently dynamic tracking and such only works for square format cameras 
            print 'dynamic tracking'
            masked_img = raw[mask_0_lo:mask_0_hi, mask_1_lo:mask_1_hi]
            masked_background = self.background[mask_0_lo:mask_0_hi, mask_1_lo:mask_1_hi]
        else:
            masked_img = raw
            masked_background = self.background
        
            
        absdiff = nim.absdiff(masked_img, masked_background)
        if self.dynamic_tracking_mask is False and self.tracking_mask is not None:
            absdiff *= tracking_mask
            
        threshold = max( 10, absdiff.max() - THRESHRANGE )
        print 'dynamic threshold: ', threshold 
        
        diffthresh = nim.threshold(absdiff, threshold)
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
                return uframe, None
            else:
                frame = Frame(raw, absdiff, diffthresh)
                return uframe, frame
            
        nblobs = blobs.max()
        if nblobs > 1:
            blobs = nim.find_biggest_blob(blobs)
            
        center = nim.center_of_blob(blobs)
        
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
    
    
