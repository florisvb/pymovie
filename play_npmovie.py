import sys


from pyglet.gl import *
from pyglet import window
from pyglet import image
from pygarrayimage.arrayimage import ArrayInterfaceImage
import time
import numpy
import numpy as np
import numpyimgproc as nim
import pickle
from optparse import OptionParser

def play_npmovie(npmovie, delay=0, magnify=1, mode='uimg'):

    w = window.Window(visible=False, resizable=True)
    
    
    arr = numpy.zeros([60,60], dtype=numpy.uint8)
    aii = ArrayInterfaceImage(arr)
    img = aii.texture
    w.width = img.width
    w.height = img.height
    w.set_visible()

    checks = image.create(32, 32, image.CheckerImagePattern())
    background = image.TileableTexture.create_for_image(checks)

    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    for uframe in npmovie.uframes:
        if uframe.uimg is not None:
            try:
                if mode == 'uimg':
                    arr = uframe.uimg
                if mode == 'wingimg2':
                    arr = uframe.wingimg2
                if mode == 'wingimg':
                    arr = uframe.wingimg
                if mode == 'flysegs':
                    arr = uframe.flysegs
                if mode == 'diffthresh':
                    arr = uframe.diffthresh
                if mode == 'full':
                    arr = npmovie.background
                    arr[uframe.indices[0]:uframe.indices[1], uframe.indices[2]:uframe.indices[3]] = uframe.uimg
            except:
                arr = None
                
            if arr is not None:
                if arr.dtype =='bool':
                    'converting'
                    arr = np.array(arr, dtype=np.uint8)*255
                    
                if magnify != 1:
                    newshape = [arr.shape[0]*magnify, arr.shape[1]*magnify]
                    arr = nim.rebin(arr, newshape)
                try:
                    aii.view_new_array(arr)
                except: # size changed
                    w.width = arr.shape[1]
                    w.height = arr.shape[0]
                    aii = ArrayInterfaceImage(arr)
        img = aii.texture
        # add some overlays:
        
        w.dispatch_events()
        
        background.blit_tiled(0, 0, 0, w.width, w.height)
        img.blit(0, 0, 0)
        w.flip()
        pyglet.graphics.draw(2, pyglet.gl.GL_LINES,('v2i', (10, 15, 30, 35)))
        time.sleep(delay) # slow down the playback
        
    w.close()

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--file", type="str", dest="file", default=None,
                        help="filename of npmovie")
    parser.add_option("--mode", type="str", dest="mode", default='uimg',
                        help="movie mode, ie. uimg, diffthresh, full, etc.")
    parser.add_option("--delay", type="int", dest="delay", default=0,
                        help="slow down the playback by x seconds")
    parser.add_option("--magnify", type="int", dest="magnify", default=1,
                        help="magnify the image by x amount (slows playback too)")
    (options, args) = parser.parse_args()
    
    if options.file is not None:
        fd = open(filename, 'r')
        npmovie = pickle.load(fd)
        fd.close()
    
    play_npmovie(npmovie, delay=options.delay, magnify=options.magnify, mode=options.mode)



    

