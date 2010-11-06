import pygmovie
import numpyimgproc as nim
import numpy as np
from scipy.ndimage.measurements import center_of_mass

THRESH_HOLD = 10
ROI_RADIUS = 20

###############

class Object:
    def __init__(self, pos, angle):
        self.pos = pos
        self.angle = angle


    

def calc_object_motion(movie):
    for i, frame in enumerate(movie.frames):

        print i

        try:
            movie.frames[i].obj.vel = (movie.frames[i+1].obj.pos - movie.frames[i-1].obj.pos) / (movie.fps*2)
            movie.frames[i].obj.speed = np.linalg.norm( movie.frames[i].obj.vel )
        except:
            movie.frames[i].obj.vel = np.array([0,0])
            movie.frames[i].obj.speed = 0.
    return movie
        
        
def calc_object_pos(movie):

    # define background, do background subtraction and thresholding
    background = nim.darken(movie.frames[0].raw, movie.frames[-1].raw)
    for i, frame in enumerate(movie.frames):

        print 'processing frames... ', i

        movie.frames[i].absdiff = nim.absdiff(frame.raw, background)
        movie.frames[i].diffthresh = nim.threshold(movie.frames[i].absdiff, THRESH_HOLD)
        
        ## TODO: right now this assumes the thresholding is perfect...!
        # find center of fly:
        center = center_of_mass(movie.frames[i].diffthresh)
        
        if 1:
            
            if center[0] == np.nan:
                center = (0,0)
            
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
            uimg = movie.frames[i].diffthresh[limlo_x:limhi_x, limlo_y:limhi_y]
            movie.frames[i].uimg = uimg
            movie.frames[i].uimg_center = center
            
            el = nim.fit_ellipse(uimg)
            center_relative_to_roi,xr,yr,theta_to_xaxis = el
            
            tmp = np.array([float(ROI_RADIUS)/2, float(ROI_RADIUS)/2])
            obj_pos = np.array(center) + np.array(center_relative_to_roi) - tmp
            
            if xr>=yr:
                angle = theta_to_xaxis
            else:
                angle = theta_to_xaxis+np.pi/2.
            
            movie.frames[i].obj = Object(obj_pos, angle)
            
            print xr,yr

    return movie        
        
if __name__ == '__main__':

    

    filename =  '/home/floris/data/SA1/C001H001S0006.avi' 
    movie = pygmovie.Movie()
    movie.fps = 5000.
    movie.load_source(filename)

    # make a mask
    mask = np.ones([movie.height, movie.width])
    mask[0:movie.height-movie.width, :] = 0
    movie.mask = mask

    # load frames
    movie.load_frames([0,20])
    movie = calc_object_pos(movie)
    movie = calc_object_motion(movie)
    
