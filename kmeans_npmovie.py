import numpy as np
import scipy.cluster.vq as cluster
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib
import sa1_movie_plots as smp
import matplotlib.patches as patches

import copy
import mdp

def get_features_for_kmeans(npmovie, frames='sync', obs='yaw,speed,slip'):
    
    frames2d = npmovie.sync2d3d.frames2d
    frames3d = npmovie.sync2d3d.frames3d

    rawyaw = npmovie.flycoord.worldangle[frames2d]
    Q = 0.1*np.eye(2) # process noise
    Q[1,1] = .01
    R = .1*np.eye(1) # observation noise
    yaw, errs = smp.smooth(rawyaw, Q,R)
    npmovie.sync2d3d.smoothyaw = yaw
    
    
    npmovie.sync2d3d.flydra_flycoord_velocities = np.zeros_like(npmovie.sync2d3d.flydra_coordinates_longaxis)
    for i in range(npmovie.sync2d3d.flydra_coordinates_longaxis.shape[0]):
        npmovie.sync2d3d.flydra_flycoord_velocities[i,0] = np.dot(npmovie.trajec.velocities[npmovie.sync2d3d.frames3d][i,0:2], npmovie.sync2d3d.flydra_coordinates_longaxis[i])
        
        npmovie.sync2d3d.flydra_flycoord_velocities[i,1] = np.sqrt( np.linalg.norm(npmovie.trajec.velocities[npmovie.sync2d3d.frames3d][i,0:2])**2 - npmovie.sync2d3d.flydra_flycoord_velocities[i,0]**2 )
    
    if obs == 'yaw,speed,slip':
        new_obs = np.nan_to_num(np.array([np.abs(yaw[:,1]), npmovie.trajec.speed[frames3d], np.abs(npmovie.flycoord.slipangle[frames2d][:,0]) ]).T)
    if obs == 'all':
        obs_list = [np.abs(yaw[:,1]), npmovie.sync2d3d.flydra_flycoord_velocities[:,0], npmovie.sync2d3d.flydra_flycoord_velocities[:,1]]
        new_obs = np.nan_to_num(np.array(obs_list).T)
        #np.abs(npmovie.flycoord.slipangle[frames2d][:,0]), npmovie.trajec.speed[frames3d]
    
    return new_obs

def find_cluster_means(movie_dataset, plot=True):
    nclusters = 5
    obs = None
    for key,npmovie in movie_dataset.items():
        if npmovie.trajec is not None:
            #print 'processing: ', key
            new_obs = get_features_for_kmeans(npmovie)
                
            if obs is None:
                obs = new_obs
            else:
                obs = np.vstack( (obs, new_obs) )
    
    nfeatures = obs.shape[1]
    obs_whitened = cluster.whiten(obs)
    feature_factors = (obs/obs_whitened)[0,:] # standard dev of the features
    cluster_means, idx = cluster.kmeans2(obs_whitened, nclusters)
    print 'idx: ', idx[0:10]
        
    cluster_errs = np.zeros( [nclusters, len(idx), nfeatures] )
    cluster_vals = np.zeros( [nclusters, len(idx), nfeatures] )
    for i, clus in enumerate(idx):
        dat = obs_whitened[i,:]
        err = dat - cluster_means[clus,:]
        cluster_errs[clus, i, :] = err
        cluster_vals[clus, i, :] = dat
    
    cluster_std = np.std(cluster_vals, axis=1)
    cluster_max = np.max( cluster_errs, axis=1)
    cluster_min = np.min( cluster_errs, axis=1)
    
    
    for key,npmovie in movie_dataset.items():
        if npmovie.trajec is not None:
            npmovie.cluster_means = cluster_means
            npmovie.feature_factors = feature_factors
            npmovie.cluster_std = cluster_std
            npmovie.cluster_max = cluster_max
            npmovie.cluster_min = cluster_min
            
    if plot:
        for key,npmovie in movie_dataset.items():
            if npmovie.trajec is not None:
                break
        plot_cluster_means(npmovie)
        
    return cluster_means, feature_factors, cluster_std
    
        
def plot_cluster_colors(npmovie, cluster_means, feature_factors, plot=True):

    obs = get_features_for_kmeans(npmovie)
    obs_whitened = obs / feature_factors
    
    cluster_assignments = np.zeros(obs.shape[0])
    
    for o in range(obs_whitened.shape[0]):
        data = obs_whitened[o,:]
        errs = np.ones(cluster_means.shape[0])*100
        for r in range(cluster_means.shape[0]):
            cluster_mean = cluster_means[r,:]
            errs[r] = np.linalg.norm(data - cluster_mean)        
        cluster_assignments[o] = np.argmin(errs)
    
    npmovie.cluster_assignments = cluster_assignments
    
    if plot:
        frames2d = npmovie.sync2d3d.frames2d
        colormap = plt.get_cmap('jet')
        colors = npmovie.cluster_assignments / float(cluster_means.shape[0])
        plt.scatter(npmovie.kalmanobj.positions[frames2d, 0], npmovie.kalmanobj.positions[frames2d, 1], 20, c=colors, cmap=colormap, edgecolors='none')
        
        
def get_cluster_for_data(npmovie, obs):
    cluster_means = npmovie.cluster_means
    feature_factors = npmovie.feature_factors
    obs_whitened = obs / feature_factors 

    errs = np.ones(cluster_means.shape[0])*100
    for r in range(cluster_means.shape[0]):
        cluster_mean = cluster_means[r,:]
        errs[r] = np.linalg.norm(obs_whitened - cluster_mean)        
    cluster = np.argmin(errs)
    return cluster
    
    
        
def plot_cluster_means(npmovie, figure=None):    
    cluster_means = npmovie.cluster_means
    cluster_min = npmovie.cluster_min
    cluster_max = npmovie.cluster_max
    cluster_std = npmovie.cluster_std
    
    nclusters = cluster_means.shape[0]
    nfeatures = cluster_means.shape[1]
    colormap = plt.get_cmap('jet')
    
    fig = pylab.figure(figure)
    theta_labels = ['yaw rate', 'speed', 'slip']
    
    for r in range(nclusters):
    
        rmax = np.max(np.abs(cluster_means))*1.2
        ax = fig.add_subplot(2,np.ceil(nclusters/2.),r+1,polar=True)
        ax.set_frame_on(False)
        
        # set cluster color
        backgroundcolor=np.array(list(colormap(r/float(nclusters))))
        theta = np.linspace(0, np.pi*2, 100)
        rad = np.ones_like(theta)*rmax
        pylab.polar( theta, rad, color=backgroundcolor, linewidth=6)
        
        theta_arr = []
        for c in range(nfeatures):
            
            theta = c / float(nfeatures)*np.pi
            theta_arr.append(theta)
            rad = cluster_means[r,c]
            
            print r,c,rad
            
            color = colormap( c / float(nfeatures))
            pylab.polar( theta, rad, 'o', color=color, markersize=7)
            
            std_line = np.linspace(rad+cluster_min[r,c], rad+cluster_max[r,c], 4, endpoint=True)
            theta_line = theta*np.ones_like(std_line)
            pylab.polar( theta_line, std_line, color=color, linewidth=.75 )
            
            std_line = np.linspace(rad-cluster_std[r,c], rad+cluster_std[r,c], 4, endpoint=True)
            theta_line = theta*np.ones_like(std_line)
            pylab.polar( theta_line, std_line, color=color, linewidth=3, )
            
        ax.set_thetagrids(np.array(theta_arr)*180/np.pi, labels=theta_labels)
        ax.set_rgrids([rmax], alpha = 0)
        ax.set_rmax(rmax)
        

def pca(data):

    nfeatures = data.shape[1]
    
    pcaanalysis = []
    
    for i in range(1,nfeatures+1):
        pcanode = mdp.nodes.PCANode()
        pcanode.output_dim = i
        pcanode.train(data)
        pcanode.stop_training()
        pcaanalysis.append( copy.copy(pcanode) )
    return pcaanalysis
        
        
        
        
        
        
        
        
        
        
        
        
        
        
