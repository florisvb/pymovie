import pickle
import sa1_movie_plots as smp
import numpy as np

def save(data, filename):
    fname = (filename)  
    fd = open( fname, mode='w' )
    pickle.dump(data, fd)
    fd.close()
    return 1

def write_to_fsee_old(npmovie, name='fsee_data'):
    
    def rotz(theta):
        ''' Returns a 3x3 rotation matrix corresponding to rotation around the *z* axis. '''
        return np.array([
            [ np.cos(theta), -np.sin(theta), 0],
            [ np.sin(theta), np.cos(theta), 0],
            [0, 0, 1]])
        
    frames = smp.get_active_frames(npmovie)
    positions = []
    attitudes = []
    states = {}
    
    flydra_coordinates_longaxis = smp.sa1_to_flydra_img_transform(npmovie.kalmanobj.long_axis[frames], fix_sign=False)*-1
    sa1_theta = np.arctan2(flydra_coordinates_longaxis[:,1], flydra_coordinates_longaxis[:,0])
    #sa1_theta = np.arctan2(npmovie.kalmanobj.long_axis[:,0], npmovie.kalmanobj.long_axis[:,1])[frames]
    sa1_time = npmovie.timestamps[frames]
    flydra_theta, flydra_time = smp.interp_sa1_data_to_flydra(npmovie, sa1_theta, sa1_time, return_time=True)
    
    for i in range(len(flydra_theta)):
        flydra_frame = i + npmovie.sa1_start_index
        positions.append( npmovie.trajec.positions[flydra_frame,:] )
        print i, flydra_frame, npmovie.trajec.positions[flydra_frame,:]
        theta = flydra_theta[i]
        attitudes.append( rotz(theta) )
        
    
    states.setdefault('positions', positions)
    states.setdefault('attitudes', attitudes)
    
    save(states, name)
    
    return 1    
    
def write_to_fsee(npmovie, name='fsee_data'):
    
    def rotz(theta):
        ''' Returns a 3x3 rotation matrix corresponding to rotation around the *z* axis. '''
        return np.array([   [ np.cos(theta), -np.sin(theta), 0],
                            [ np.sin(theta), np.cos(theta), 0],
                            [0, 0, 1]])
        
    frames = smp.get_active_frames(npmovie)
    positions = []
    attitudes = []
    states = {}
    
    for i in range(len(npmovie.sync2d3d.pos2d)):
        positions.append( npmovie.sync2d3d.pos3d_est[i,:] )
        print npmovie.sync2d3d.yaw[i]
        attitudes.append( rotz(npmovie.sync2d3d.yaw[i]) )
    
    states.setdefault('positions', positions)
    states.setdefault('attitudes', attitudes)
    
    save(states, name)
    
    return 1 
