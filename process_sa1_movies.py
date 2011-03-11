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
import flydra_floris_analysis as ffa
import classification as clas
import copy


def test(n):
    
    o, m = sa1a.sa1_analysis()
    print 'processing: ', m.keys()[n]
    dataset = pm.load('/home/floris/data/sa1_movie_data/h5_files/reduced_dataset')
    npmovie = process_movie(m.keys()[n], m, dataset)
    
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
        
        npmovie.id = MOVIE_ID
        npmovie.behavior = movie_info[MOVIE_ID]['Behavior']
        npmovie.path = movie_info[MOVIE_ID]['Path']
        npmovie.posttype = movie_info[MOVIE_ID]['PostType']
        npmovie.extras = movie_info[MOVIE_ID]['Extras']
        npmovie.trigger_stamp = movie_info[MOVIE_ID]['Trigger Stamp']
        npmovie.timestamps = movie_info[MOVIE_ID]['Timestamps']

        # make a mask
        mask = np.ones([movie.height, movie.width])
        mask[0:movie.height-movie.width, :] = 0
        movie.mask = mask

        # make a tracking mask -- 0 where no tracking desired
        tracking_mask = np.abs(nim.plot_circle(1024,1024,[512,512], 150)-1)
        npmovie.tracking_mask = tracking_mask

        npmovie.calc_background(movie)
        npmovie.load_frames(movie, timerange=None)
        npmovie.fps = 5000.
        npmovie.obj = pm.Object(npmovie)
        pm.auto_adjust_uimg(npmovie)
        pm.segment_fly(npmovie)
        
        frames = smp.get_active_frames(npmovie)
        if len(frames) < 5: 
            del(movie)
            return npmovie
        
        pm.calc_obj_pos(npmovie)
        pm.calc_obj_motion(npmovie)
        pm.smooth_legs(npmovie)
        pm.calc_polarpos(npmovie)
        pm.smooth_axis_ratio(npmovie)
        
        trajec, time_err = sa1a.get_flydra_trajectory(npmovie, dataset)
        print 'trajectory in process_movie: ', trajec
        npmovie.trajec = copy.copy(trajec)
            
        clas.calc_fly_coordinates(npmovie)
        smp.calc_frame_of_landing(npmovie)
        
        if trajec is not None:
            smp.sync_2d_3d_data(npmovie)
            if len(npmovie.sync2d3d.frames3d) > 2:
                smp.calc_body_orientations(npmovie)
                smp.calc_smoothyaw(npmovie, plot=False)
        else:
            npmovie.sync2d3d = None
        
        npmovie.cluster = None
        #strobe_img = pm.strobe_image(movie)
        #movie_info[MOVIE_ID].setdefault('StrobeImg', strobe_img)

        #npmovie_filename = '/home/floris/data/windtunnel/SA1/sa1_npmovie_' + MOVIE_ID
        #pm.save(npmovie, npmovie_filename)
        
        del(movie)
        return npmovie
        
def get_relevant_obj_ids_from_dataset(dataset, obj_id_list):
    reduced_dataset = ffa.Dataset(like=dataset)
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
                        reduced_dataset.trajecs.setdefault(copy.copy(k), copy.copy(trajectory))
                else:
                    continue
            
    return reduced_dataset                
                
                
        
def batch_process_harddrive(path='/media/SA1_videos/sa1_movies_3/', movie_info=None, dataset=None):

    movie_dataset_filename = None

    if movie_info is None:
        obj_id_list, movie_info = sa1a.sa1_analysis()
    
    for MOVIE_ID in movie_info.keys():
        #if movie_info[MOVIE_ID]['Path'] == path:
        print 'loading id: ', MOVIE_ID

        npmovie = process_movie(MOVIE_ID, movie_info, dataset=dataset)
        print 'trajectory: ', npmovie.trajec
        mnpmovie = pm.MiniNPM(npmovie)
        print 'mini trajectory: ', mnpmovie.trajec
        del(npmovie)
        
        if movie_dataset_filename is not None:
            movie_dataset = pm.load(movie_dataset_filename)
        else:
            movie_dataset = {}
            movie_dataset_filename = 'movie_dataset_3'
        movie_dataset.setdefault(MOVIE_ID, mnpmovie)
        pm.save(movie_dataset, movie_dataset_filename)
        del(movie_dataset)
              
def get_trajecs_for_movie_dataset(movie_dataset, movie_info, dataset):  
    
    for k, npmovie in movie_dataset.items(): 
        print npmovie.id
        npmovie.objid = movie_info[npmovie.id]['Obj_ID']

        print k
        trajec, time_err = sa1a.get_flydra_trajectory(npmovie, dataset)
        print trajec
        npmovie.trajec = copy.copy(trajec)
        
    return movie_dataset
    
        
if __name__ == '__main__':        

    dataset = ap.initialize_dataset('/home/floris/data/windtunnel/SA1/dataset_all')
    batch_process_harddrive(path='/media/SA1_movies_3/sa1_movies/', movie_dataset_filename = None, dataset=dataset)


