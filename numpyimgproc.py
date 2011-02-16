import numpy as np
import matplotlib.pyplot as plt

import scipy.optimize
from scipy.ndimage.measurements import center_of_mass
import scipy.ndimage as ndimage

import numpy as n
import scipy.interpolate
from scipy.ndimage.morphology import binary_erosion

inf = np.inf

def in_range(val, rang):
    if val > rang[0] and val < rang[1]:
        return True
    else:
        return False

def compare(a, b=None, method=None):

    if type(a) is list:
        result = a[0]
        for i in range([1,len(a)]):
            if method == 'lighten':
                result = np.maximum(result,a[i])
            elif method == 'darken':
                result = np.minimum(result,a[i])
        return result
    elif b is not None:
        if method is 'lighten':
            result = np.maximum(a,b)
        elif method is 'darken':
            result = np.minimum(a,b)
        return result
    else:
        ValueError('please enter a valid a,b pair')
        
def darken(a, b=None):
    result = compare(a,b,method='darken')
    return result
def lighten(a, b=None):
    result = compare(a,b,method='lighten')
    return result
    
def threshold(img, threshold_lo, threshold_hi=255):
    threshed_lo = img>threshold_lo
    threshed_hi = img<threshold_hi
    threshed = threshed_lo*threshed_hi
    threshed *= 255
    return threshed
    
def absdiff(a,b):
    img = np.array(a, dtype=float)
    bkg = np.array(b, dtype=float)
    diff = img-bkg
    absdiff = np.abs(diff)
    return absdiff
    
def auto_adjust_levels(img):
      
    img2 = (img-img.min())
    if img2.max() > 0:
        img3 = img2*int(255/float(img2.max()))
    else:
        img3 = img2
        
    return img3
    
def find_blobs(img, sizerange=[0,inf], aslist=False):
    blobs, nblobs = ndimage.label(img)
    blob_list = []
    if nblobs < 1:
        if aslist is False:
            return np.zeros_like(img)
        else:
            return [np.zeros_like(img)]
    #print 'n blobs: ', nblobs
    # erode filter
    n_filtered = 0
    for n in range(1,nblobs+1):
        blob_size = (blobs==n).sum()
        if not in_range(blob_size, sizerange):
            blobs[blobs==n] = 0
            nblobs -= 1
        else:
            if aslist:
                b = np.array(blobs==n, dtype=np.uint8)
                blob_list.append(b)
            else:
                n_filtered += 1
                blobs[blobs==n] = n_filtered
    
    if aslist is False:
        if nblobs < 1:
            return np.zeros_like(img)
        else:
            blobs = np.array( (blobs>0)*255, dtype=np.uint8)
            return blobs
    else:
        if len(blob_list) < 1:
            blob_list = [np.zeros_like(img)]
        return blob_list
        
def find_biggest_blob(img):
    blobs, nblobs = ndimage.label(img)
    
    if nblobs < 1:
        return None
    
    if nblobs == 1:
        return blobs

    biggest_blob_size = 0
    biggest_blob = None
    for n in range(1,nblobs+1):
        blob_size = (blobs==n).sum()
        if blob_size > biggest_blob_size:
            biggest_blob = blobs==n
            biggest_blob_size = blob_size
    biggest_blob = np.array(biggest_blob, dtype=np.uint8)
    return biggest_blob 
    
def find_2_biggest_blobs(img):
    blobs, nblobs = ndimage.label(img)
    b1 = find_biggest_blob(blobs)
    blobs[b1] = 0
    b2 = find_biggest_blob(blobs)
    return b1, b2
    
def find_n_biggest_blobs(img, n=1):
    blobs, nblobs = ndimage.label(img)
    
    blob_list = []
    for i in range(n):
        b = find_biggest_blob(blobs)
        if b == None:
            break
        blob_list.append(b)
        blobs[b>0] = 0
    return blob_list

def find_smallest_blob(img):
    blobs, nblobs = ndimage.label(img)
    
    if nblobs < 1:
        return None
    
    if nblobs == 1:
        return blobs

    smallest_blob_size = 10000
    smallest_blob = None
    for n in range(1,nblobs+1):
        blob_size = (blobs==n).sum()
        if blob_size < smallest_blob_size:
            smallest_blob = blobs==n
            smallest_blob_size = blob_size
            
    return smallest_blob 


def plot_circle(xmesh, ymesh, center, r):
    center = list(center)
    
    def in_circle(x,y):
        R = np.sqrt((X-center[1])**2 + (Y-center[0])**2)
        Z = R<r
        return Z

    x = np.arange(0, xmesh, 1)
    y = np.arange(0, ymesh, 1)
    X,Y = np.meshgrid(x, y)

    Z = in_circle(X, Y)
    
    return Z
    
def fit_circle(img):

    center = center_of_mass(img)
    r = 10
    
    x = np.arange(0, img.shape[0], 1)
    y = np.arange(0, img.shape[1], 1)
    X,Y = np.meshgrid(x, y)    

    def circle(x,y, center, r):
        R = np.sqrt((X-center[1])**2 + (Y-center[0])**2)
        Z = R<r
        Z *= 255
        return Z

    def fminfunc(vals):
        #center = [vals[0], vals[1]]
        r = vals[0]
        Z = circle(X, Y, center, r)
        err = np.abs(img-Z)
        errsum = np.sum(err)
        #print errsum
        return errsum
    
    vals = [r]
    tmp = scipy.optimize.fmin( fminfunc, vals, full_output = 1, disp=0)
    vals = tmp[0]
    #center = [vals[0], vals[1]]
    r = vals[0]
    Z = circle(X, Y, center, r)
    
    return center, r, Z
    
def ellipse(x,y,center,a,b,angle):

    a = np.abs(a)
    b = np.abs(b)
    angle = np.abs(angle)    
    
    t = np.tan(angle)
    xaxis = [1,t]
    xaxis /= np.linalg.norm(xaxis)
    Xa, Ya = xaxis
    Xc = center[1]
    Yc = center[0]
    
    xe = Xa*(x - Xc) + Ya*(y - Yc)
    ye = -1*Ya*(x  - Xc) + Xa*(y - Yc)
    
    R = (xe)**2/a**2 + (ye)**2/b**2
    Z = R<1
    Z *= 255.
    return Z


def fit_ellipse(img, full_output=False):

    img = np.array(img, dtype=float)

    center = center_of_mass(img)
    r = 10.
    
    x = np.arange(0., float(img.shape[1]), 1.)
    y = np.arange(0., float(img.shape[0]), 1.)
    X,Y = np.meshgrid(x, y)    

    def fminfunc(vals):
        #center = [vals[0], vals[1]]
        a,b,angle = vals
        Z = ellipse(X, Y, center, a, b, angle)
        err = np.abs(img-Z)
        errsum = np.sum(err)
        #print errsum
        return errsum
    
    a = 20
    b = 20
    angle = np.pi/4
    vals = [a,b,angle]
    tmp = scipy.optimize.fmin( fminfunc, vals, full_output = 1, disp=0)
    vals = tmp[0]
    #center = [vals[0], vals[1]]
    a,b,angle = vals
    
     
    if full_output is False:
        return center, a, b, angle
    else:
        Z = ellipse(X, Y, center, a, b, angle)
        return center, a, b, angle, Z
        
def center_of_blob(img):
    if type(img) is list:
        centers = []
        for blob in img:
            center = np.array([center_of_mass(blob)[i] for i in range(2)])
            centers.append(center)
        return centers
    else:
        center = np.array([center_of_mass(img)[i] for i in range(2)])
        return center
    
def extract_blob(img, radius=None):
    center = center_of_blob(img)
    pts = np.transpose(np.nonzero(img)).T
    if radius is None: # find minumum bounding box
        x_lo = pts[0,:].min()
        x_hi = pts[0,:].max()
        y_lo = pts[1,:].min()
        y_hi = pts[1,:].max()
        radius = int(max( center[0]-x_lo, x_hi-center[0], center[1]-y_lo, y_hi-center[1] ))
        print center, radius
        print x_lo, x_hi, y_lo, y_hi
        
    uimg = img[center[0]-radius:center[0]+radius+1, center[1]-radius:center[1]+radius+1]
    return center, uimg
        
def fit_ellipse_cov(img, erode=True, recenter=False):
    # Pattern. Recogn. 20, Sept. 1998, pp. 31-40
    # J. Prakash, and K. Rajesh
    # Human Face Detection and Segmentation using Eigenvalues of Covariance Matrix, Hough Transform and Raster Scan Algorithms

    #eroded_img = binary_erosion(img)
    #boundary = img-eroded_img
    
    if img is not None:
    
        if erode:
            try:
                img = binary_erosion(img)
            except:
                pass
                
        if recenter:
            center = center_of_blob(img)
        else:
            center = np.array([0,0])

        if 1:
            pts = np.transpose(np.nonzero(img))
            for pt in pts:
                pt -= center
            pts = (np.transpose(np.nonzero(img))).T
            cov = np.cov(pts)
            cov = np.nan_to_num(cov)
            e,v = np.linalg.eig(cov)
            ratio = max(e) / min(e)
            
            i = np.argmax(e)
            long_axis = v[:,i]
        if recenter is False:
            return long_axis, ratio
        else:
            return center, long_axis, ratio
            
    else:
        return [0,0],0
        
############################################################################


def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = n.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None

def place_images(image_size, images, image_placements):
    # this function is used to put separate images together on one image plane in arbitrary placements - for making composite movies
    # currently assumes grayscale!
    master_image = np.zeros(image_size)
    for i, image in enumerate(images):
        if len(image.shape) > 2:
            image = image[:,:,0]
        if image.max() > 0:
            image = np.array(image, dtype=float)
            image /= image.max()
        origin = image_placements[i]
        master_image[origin[0]:origin[0]+image.shape[0], origin[1]:origin[1]+image.shape[1]] = image
    master_image_rgb = np.zeros([master_image.shape[0], master_image.shape[1], 3])
    for i in range(3):
        master_image_rgb[:,:,i] = master_image
    return master_image_rgb    
    
def rotate_image(img, rot):
    
    imgrot = np.zeros_like(img)
    
    for r in range(img.shape[0]):
        for w in range(img.shape[1]):
            ptrot = np.dot(rot, np.array([r,w]))
            rrot = ptrot[0]
            wrot = ptrot[1]
            if rrot < 0:
                rrot += img.shape[0]
            if wrot < 0:
                wrot += img.shape[1]
            imgrot[rrot, wrot] = img[r,w]
    return imgrot
                
    
##############################################################################
def rebin( a, newshape ):
        '''Rebin an array to a new shape.
        '''
        assert len(a.shape) == len(newshape)

        slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
        coordinates = scipy.mgrid[slices]
        indices = coordinates.astype('i')   #choose the biggest smaller integer index
        return a[tuple(indices)]

