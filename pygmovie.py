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
import cPickle as pickle
import copy

import numpy
import adskalman.adskalman as adskalman

THRESHHOLD = 20
WINGTHRESHOLD = 40
ROI_RADIUS = 30
############################################################################################

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
    
    for uframe in npmovie.uframes:
        if uframe is not None:
            uimg = uframe.uimg
            uimg_no_body = uframe.uimg + uframe.diffthresh*255
            uimg_no_body_thresh = uimg_no_body<WINGTHRESHOLD
            wingimg = binary_erosion(uimg_no_body_thresh)
            uframe.wingimg = wingimg
            
            
            
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
    



############################################################################################
class Object:
    def __init__(self, npmovie):
        self.positions = np.zeros([len(npmovie.uframes), 2])
        self.velocities = np.zeros([len(npmovie.uframes), 2])
        self.long_axis = np.zeros([len(npmovie.uframes), 2])
        

def calc_obj_motion(npmovie, fps):
    npmovie.fps = fps
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
    

    y = np.vstack( (npmovie.obj.positions[:,0], npmovie.obj.positions[:,1]) ).T
    
    initx = numpy.array([npmovie.obj.positions[0,0], npmovie.obj.positions[0,1],0,0],dtype=numpy.float)
    initV = 0*numpy.eye(ss)

    xsmooth,Vsmooth = adskalman.kalman_smoother(y,F,H,Q,R,initx,initV)

    npmovie.kalmanobj.positions = xsmooth
    npmovie.kalmanobj.velocities = Vsmooth
    
    return npmovie


    
        
        
def calc_obj_pos(npmovie):

    for i in range(len(npmovie.uframes)):
        uframe = npmovie.uframes[i]
        if uframe is not None:
            
            npmovie.obj.positions[i,:] = uframe.center
            npmovie.obj.long_axis[i,:] = nim.fit_ellipse_cov(uframe.diffthresh)
                
    return npmovie
################################################################################################################

class Movie:
    def __init__(self, filename=None):
        self.mask = None
        if filename is not None:
            self.load_source(filename)
    
    def load_source(self, filename):
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
        
        self.blob_size_range = [100,250]
        
        self.fps = None
        
    def calc_background(self, movie):
        
        first_frame = movie.get_frame_at_timestamp(0)
        last_frame = movie.get_frame_at_timestamp(movie.duration-1)
        self.background = nim.lighten(first_frame, last_frame)
        self.height, self.width = self.background.shape
        
    def background_subtraction(self, raw, save_raw=False):

        if self.dynamic_tracking_mask is True:
            print 'dynamic tracking'
            mask_center_0 = int(self.mask_center[0])
            mask_center_1 = int(self.mask_center[1])
            mask_0_lo = max(0, mask_center_0 - self.mask_radius)
            mask_0_hi = min(self.width, mask_center_0 + self.mask_radius)
            mask_1_lo = max(0, mask_center_1 - self.mask_radius)
            mask_1_hi = min(self.width, mask_center_1 + self.mask_radius)
            masked_img = raw[mask_0_lo:mask_0_hi, mask_1_lo:mask_1_hi]
            masked_background = self.background[mask_0_lo:mask_0_hi, mask_1_lo:mask_1_hi]
        else:
            masked_img = raw
            masked_background = self.background
        
            
        absdiff = nim.absdiff(masked_img, masked_background)
        diffthresh = nim.threshold(absdiff, THRESHHOLD)
            
        if self.dynamic_tracking_mask is False and self.tracking_mask is not None:
            diffthresh *= self.tracking_mask
            
        
        ## TODO: right now this assumes the thresholding is perfect...!
        # find center of fly:
        blobs = nim.find_blobs(diffthresh, self.blob_size_range)
        
        if blobs is None:
            uframe = uFrame()
            if save_raw is False:
                return uframe, None
            else:
                frame = Frame(raw, absdiff, diffthresh)
                return uframe, frame
            
        nblobs = blobs.max()
        if nblobs > 1:
            blobs = nim.find_smallest_blob(blobs)
            
        center = center_of_mass(blobs)
        
        if 1:
            limlo_x = center[0]-ROI_RADIUS
            if limlo_x < 0:
                limlo_x = 0
            limhi_x = center[0]+ROI_RADIUS
            if limhi_x > movie.height:
                limhi_x = movie.height
            limlo_y = center[1]-ROI_RADIUS
            if limlo_y < 0:
                limlo_y = 0
            limhi_y = center[1]+ROI_RADIUS
            if limhi_y > movie.width:
                limhi_y = movie.width
                
            # TODO: right now object entering or leaving frame doesnt work perfectly
            uimg = masked_img[limlo_x:limhi_x, limlo_y:limhi_y]
            uimg = copy.copy(uimg)
            uimg_absdiff = absdiff[limlo_x:limhi_x, limlo_y:limhi_y]
            uimg_absdiff = copy.copy(uimg_absdiff)
            uimg_diffthresh = blobs[limlo_x:limhi_x, limlo_y:limhi_y]
            uimg_diffthresh = copy.copy(uimg_diffthresh)
            
            if self.dynamic_tracking_mask is True:
                tmp = np.array([float(self.mask_radius), float(self.mask_radius)])
                center = np.array(self.mask_center) + np.array(center) - tmp
            self.dynamic_tracking_mask = True
            self.mask_center = center
            
            print center
            
            uframe = uFrame(center, uimg, uimg_absdiff, uimg_diffthresh) 
            uframe.blobsize = blobs.sum()
                
            if save_raw is False:
                del(raw)
                del(absdiff)
                del(diffthresh)
                del(blobs)
                return uframe, None
                
            frame = Frame(raw, absdiff, diffthresh)
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
                
            raw = movie.get_next_frame()    
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
            print 'loading... frame num: ', fr, '  ::  ', 'timestamp: ', timestamp, ' :: ', timestamp/movie.duration*100, '%'
            
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
        self.timestamp = None
    def show(self, var='uimg'):
        if var == 'uimg':
            plt.imshow(self.uimg)
        elif var == 'absdiff':
            plt.imshow(self.absdiff)
        elif var == 'diffthresh':
            plt.imshow(self.diffthresh)  
        else:
            print 'please use a valid var'
            return
        plt.show()
        return        
        
if __name__ == '__main__':

    filename =  '/media/SA1_videos/sa1_videos/20101027_C001H001S0006.avi'
    
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
    
    save(npmovie, 'sa1_movie')
    
    npmovie.obj = Object(npmovie)
    calc_obj_pos(npmovie)
    calc_obj_motion(npmovie, 5000.)
    save(npmovie, 'sa1_movie_obj')
    
    
