import sys


from pyglet.gl import *
from pyglet import window
from pyglet import image
from pygarrayimage.arrayimage import ArrayInterfaceImage
import time
import numpy
import numpyimgproc as nim

def play_npmovie(npmovie, delay=0, magnify=1):

    w = window.Window(visible=False, resizable=True)
    
    
    arr = numpy.zeros([600,600], dtype=numpy.uint8)
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
        try:
            arr = uframe.uimg
            if magnify != 1:
                newshape = [arr.shape[0]*magnify, arr.shape[1]*magnify]
                arr = nim.rebin(arr, newshape)
            try:
                aii.view_new_array(arr)
            except: # size changed
                aii = ArrayInterfaceImage(arr)
                w.width = img.width
                w.height = img.height
                w.set_visible()
        except:
            pass
        img = aii.texture
        
        w.dispatch_events()
        
        background.blit_tiled(0, 0, 0, w.width, w.height)
        img.blit(0, 0, 0)
        w.flip()
        
        time.sleep(delay) # slow down the playback
        
    w.close()

