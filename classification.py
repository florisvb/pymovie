import numpy as np
import copy

class Object:
    def __init__(self):
        pass

def rotate_coordinates(old_coordinates, basis1, basis2):
    

    if len(old_coordinates.shape) == 1:
        b1 = basis1
        b2 = basis2
        n = np.zeros_like(old_coordinates)
        n[0] = n1 = np.dot(old_coordinates, b1)
        n[1] = np.dot(old_coordinates, b2)
        return n

    elif len(old_coordinates.shape) == 2:
    
        n = np.zeros_like(old_coordinates)
        for i in range(len(old_coordinates)):
            b1 = basis1[i]
            b2 = basis2[i]
            n[i,0] = np.dot(old_coordinates[i], b1)
            n[i,1] = np.dot(old_coordinates[i], b2)
        return n
        
def iseven(n):
    if int(n)/2.0 == int(n)/2:
        return True
    else:
        return False 
       
def isodd(n):
    if int(n)/2.0 == int(n)/2:
        return False
    else:
        return True
        
def remove_angular_rollover(A, max_change_acceptable):
    array = copy.copy(A)
    for i, val in enumerate(array):
        if i == 0:
            continue
            
        diff = array[i] - array[i-1]
        if np.abs(diff) > max_change_acceptable:
            factor = np.round(np.abs(diff)/(np.pi))  
            if iseven(factor):
                array[i] -= factor*np.pi*np.sign(diff)
                
    if len(A) == 2:
        return array[1]
    else:
        return array
    
def calc_fly_coordinates(npmovie, post_pos = [512, 512]):
    post_pos = np.array(post_pos)
    npmovie.flycoord = Object()
    
    
    npmovie.kalmanobj.short_axis = np.zeros_like(npmovie.kalmanobj.long_axis)
    npmovie.flycoord.heading = np.zeros_like(npmovie.kalmanobj.long_axis)
    npmovie.flycoord.slipangle = np.zeros_like(npmovie.kalmanobj.speed)
    npmovie.flycoord.postangle = np.zeros_like(npmovie.kalmanobj.speed)
    npmovie.flycoord.worldangle = np.zeros_like(npmovie.kalmanobj.speed)
    npmovie.flycoord.speed = npmovie.kalmanobj.speed
    npmovie.flycoord.dist_to_post = np.zeros_like(npmovie.kalmanobj.speed)
    
    for i in range(len(npmovie.kalmanobj.short_axis)):
        npmovie.kalmanobj.short_axis[i,0] = npmovie.kalmanobj.long_axis[i,1]
        npmovie.kalmanobj.short_axis[i,1] = npmovie.kalmanobj.long_axis[i,0]*-1
        
        npmovie.flycoord.worldangle[i] = np.arctan2(npmovie.kalmanobj.long_axis[i,1], npmovie.kalmanobj.long_axis[i,0]) # want angle between 0 and 360 deg
        if npmovie.flycoord.worldangle[i] < 0:
            npmovie.flycoord.worldangle[i] = np.pi*2+npmovie.flycoord.worldangle[i]
        # smooth world angle...
        
        if 0:
            eps = .5
            if i>0:
                diff = npmovie.flycoord.worldangle[i] - npmovie.flycoord.worldangle[i-1]
                if np.abs(diff) > np.pi*2-eps and np.abs(diff) < np.pi*4-eps:
                    npmovie.flycoord.worldangle[i] -=  2*np.pi*np.sign(diff)
                if np.abs(diff) > np.pi*4-eps and np.abs(diff) < np.pi*6-eps:
                    npmovie.flycoord.worldangle[i] -=  2*np.pi*np.sign(diff)
    
        if 1:
            if i>0:
                npmovie.flycoord.worldangle[i] = remove_angular_rollover( [npmovie.flycoord.worldangle[i-1], npmovie.flycoord.worldangle[i]], .5 )
            
        npmovie.flycoord.heading[i] = npmovie.kalmanobj.velocities[i] / np.linalg.norm( npmovie.kalmanobj.velocities[i] )
        cosslipangle = np.dot(npmovie.flycoord.heading[i], npmovie.kalmanobj.long_axis[i]) / (np.linalg.norm(npmovie.flycoord.heading[i])*np.linalg.norm(npmovie.kalmanobj.long_axis[i])) 
        npmovie.flycoord.slipangle[i] = np.arccos(cosslipangle)
        vec_to_post = post_pos - npmovie.kalmanobj.positions[i]
        npmovie.flycoord.dist_to_post[i] = np.linalg.norm(vec_to_post)
        cospostangle = np.dot(vec_to_post, npmovie.kalmanobj.long_axis[i]) / ( np.linalg.norm(vec_to_post)*np.linalg.norm(npmovie.kalmanobj.long_axis[i]) )
        npmovie.flycoord.postangle[i] = np.arccos(cospostangle)
            
    npmovie.flycoord.velocities = rotate_coordinates(npmovie.kalmanobj.velocities, npmovie.kalmanobj.long_axis, npmovie.kalmanobj.short_axis)
    

def calc_fly_coordinates_dataset(movie_dataset):
    
    for key,mnpmovie in movie_dataset.items():
        calc_fly_coordinates(mnpmovie)
    








    
