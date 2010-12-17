import pickle
import pygmovie as pm
import numpyimgproc as nim
import sys
import numpy as np

sys.path.append('/home/floris/src/analysis')

import analysis_plot as ap


def process_movie(filename): 
        
        movie_filename = movie_info[filename]['Path'] + filename + '.avi'
        print 'loading: ', movie_filename
        movie = pm.Movie(movie_filename)
        npmovie = pm.npMovie()

        # make a mask
        mask = np.ones([movie.height, movie.width])
        mask[0:movie.height-movie.width, :] = 0
        movie.mask = mask

        # make a tracking mask -- 0 where no tracking desired
        tracking_mask = np.abs(nim.plot_circle(1024,1024,[512,512], 150)-1)
        npmovie.tracking_mask = tracking_mask

        npmovie.calc_background(movie)
        npmovie.load_frames(movie)
        
        npmovie.fps = 5000.
        
        npmovie.obj = pm.Object(npmovie)
        pm.calc_obj_pos(npmovie)
        pm.calc_obj_motion(npmovie)
                    
        pm.find_wings(npmovie)
        pm.split_wings_blobs(npmovie)
        pm.calc_wing_motion(npmovie)
        pm.interp_nan_values(npmovie)
        pm.calc_wing_angle(npmovie)
        pm.kalmanize_wing_angle(npmovie)
            
        movie_info[filename].setdefault('KalmanObj', npmovie.kalmanobj)
        movie_info[filename].setdefault('Obj', npmovie.obj)
        
        strobe_img = pm.strobe_image(movie)
        movie_info[filename].setdefault('StrobeImg', strobe_img)

        npmovie_filename = '/home/floris/data/windtunnel/SA1/sa1_npmovie_' + filename
        pm.save(npmovie, npmovie_filename)
        
        return npmovie, movie_info



ap.initialize_dataset('/home/floris/data/windtunnel/SA1/dataset_all')

filename = '/home/floris/data/windtunnel/SA1/movie_info'
fname = (filename)  
fd = open( fname, mode='r' )
movie_info = pickle.load(fd)
fd.close()

npmovie, movie_info = process_movie('20101101_C001H001S0019')


if 0:
    f = 0
    for filename in movie_info.keys():
        if movie_info[filename]['Behavior'] == 'landing':
            #if f >= 1:
            #    break
            #f += 1
            
            npmovie, movie_info = process_movie(filename)
        
    filename = '/home/floris/data/windtunnel/SA1/movie_info_with_trajectories'
    fname = (filename)  
    fd = open( fname, mode='w' )
    pickle.dump(movie_info, fd)
    fd.close()
