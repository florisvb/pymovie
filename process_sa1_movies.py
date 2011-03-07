import pickle
import pygmovie as pm
import numpyimgproc as nim
import sys
import numpy as np
import copy
import sa1_movie_plots as smp

sys.path.append('/home/floris/src/analysis')
import sa1_analysis as sa1a
import analysis_plot as ap

def test():
    
    o, m = sa1a.sa1_analysis()
    print 'processing: ', m.keys()[-1]
    dataset = pm.load('/home/floris/data/sa1_movie_data/h5_files/reduced_dataset')
    npmovie = process_movie(m.keys()[-1], m, dataset)
    
    return npmovie

        
def process_movie(MOVIE_ID, movie_info=None, dataset=None): 
        
        if movie_info is None:
            filename = '/home/floris/data/windtunnel/SA1/movie_info'
            fname = (filename)  
            fd = open( fname, mode='r' )
            movie_info = pickle.load(fd)
            fd.close()
        
        movie_filename = movie_info[MOVIE_ID]['Path'] + MOVIE_ID + '.avi'
        print 'loading: ', movie_filename
        movie = pm.Movie(movie_filename)
        npmovie = pm.npMovie()
        
        npmovie.objid = movie_info[MOVIE_ID]['Obj_ID']
        npmovie.epochtime = movie_info[MOVIE_ID]['EpochTime']
        npmovie.filename = MOVIE_ID

        # make a mask
        mask = np.ones([movie.height, movie.width])
        mask[0:movie.height-movie.width, :] = 0
        movie.mask = mask

        # make a tracking mask -- 0 where no tracking desired
        tracking_mask = np.abs(nim.plot_circle(1024,1024,[512,512], 150)-1)
        npmovie.tracking_mask = tracking_mask

        npmovie.calc_background(movie)
        npmovie.load_frames(movie, timerange=[.1,.15])
        
        npmovie.fps = 5000.
        
        npmovie.obj = pm.Object(npmovie)
        
        pm.auto_adjust_uimg(npmovie)
        pm.segment_fly(npmovie)
        
        pm.calc_obj_pos(npmovie)
        print npmovie.obj.positions.shape
        pm.calc_obj_motion(npmovie)
        
        
        
        # wing kinematics.. need to solve the wing axis / angle problems
        if 0:
            pm.find_wings(npmovie)
            pm.split_wings_bodyaxis(npmovie)
            pm.split_wings_blobs(npmovie)
            pm.calc_wing_motion(npmovie)
            pm.interp_nan_values_R(npmovie)
            pm.interp_nan_values_L(npmovie)
            pm.calc_wing_angle(npmovie)
            pm.kalmanize_wing_angle(npmovie)
            pm.find_legs(npmovie)
            
        pm.smooth_legs(npmovie)
        pm.calc_polarpos(npmovie)
        pm.smooth_axis_ratio(npmovie)
        trajec, time_err = sa1a.get_flydra_trajectory(npmovie, dataset)
        smp.sync_2d_3d_data(npmovie)
        #strobe_img = pm.strobe_image(movie)
        #movie_info[MOVIE_ID].setdefault('StrobeImg', strobe_img)

        #npmovie_filename = '/home/floris/data/windtunnel/SA1/sa1_npmovie_' + MOVIE_ID
        #pm.save(npmovie, npmovie_filename)
        
        del(movie)
            
        return npmovie
        
def get_relevant_obj_ids_from_dataset(dataset, obj_id_list):
    reduced_dataset = {}
    for i in range(obj_id_list.shape[0]):
        obj_id = obj_id_list[i,0]
        epochtime = obj_id_list[i,1]
        
        if obj_id is not None:
            obj_id = int(obj_id)
            
            for k, trajectory in dataset.trajecs.items():
                t = k.lstrip('1234567890')
                t = t.lstrip('_')
                o = int(t)
                if obj_id == o:
                    # check timestamp:
                    time_err = np.abs(trajectory.epoch_time[0] - epochtime)
                    print time_err
                    if time_err < 100:
                        reduced_dataset.setdefault(copy.copy(k), copy.copy(trajectory))
                else:
                    continue
            
    return reduced_dataset                
                
                
        
def batch_process_harddrive(path='/media/SA1_videos/sa1_movies/', movie_dataset_filename = None, dataset=None):
    
    filename = '/home/floris/Documents/data/movie_info'
    fname = (filename)  
    fd = open( fname, mode='r' )
    movie_info = pickle.load(fd)
    fd.close()
    print len(movie_info.keys())
    
    
    for MOVIE_ID in movie_info.keys():
        #if movie_info[MOVIE_ID]['Path'] == path:
        print 'loading id: ', MOVIE_ID

        npmovie = process_movie(MOVIE_ID, movie_info)
        mnpmovie = pm.MiniNPM(npmovie)
        del(npmovie)
        
        mnpmovie.behavior = movie_info[MOVIE_ID]['Behavior']
        mnpmovie.path = movie_info[MOVIE_ID]['Path']
        mnpmovie.posttype = movie_info[MOVIE_ID]['PostType']
        mnpmovie.extras = movie_info[MOVIE_ID]['Extras']
        mnpmovie.id = MOVIE_ID
        
        if movie_dataset_filename is not None:
            movie_dataset = pm.load(movie_dataset_filename)
        else:
            movie_dataset = {}
            movie_dataset_filename = 'movie_dataset'
        movie_dataset.setdefault(MOVIE_ID, mnpmovie)
        pm.save(movie_dataset, movie_dataset_filename)
        del(movie_dataset)
              
        
if __name__ == '__main__':        

    dataset = ap.initialize_dataset('/home/floris/data/windtunnel/SA1/dataset_all')
    batch_process_harddrive(path='/media/SA1_movies_3/sa1_movies/', movie_dataset_filename = None, dataset=dataset)


