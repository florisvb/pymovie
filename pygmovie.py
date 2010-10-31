import matplotlib.pyplot as plt
import numpy
from pyglet import media
import numpy as np
import pyglet
    
class Movie:
    def __init__(self):
        self.frames = []
    def load_using_pyglet(self, filename, nframes=[0,np.inf]):
    
        try:
            player = media.Player()
            source = media.load(filename)
            player.queue(source)
            player.play()
        except:
            ValueError('failed to open movie!')
            
        # load all requested frames
        fr = nframes[0]
        while fr < nframes[1]:
            player.dispatch_events()
            imdata=player.texture.get_image_data()
            a = numpy.frombuffer(player.texture.get_image_data().data, numpy.uint8)
            a.shape = (imdata.width, imdata.height, 4)
            a = a[:,:,0:3] #drop alpha channel

            mono = (a[:,:,0]+a[:,:,1]+a[:,:,2])/3 #convert to monochrome
            frame = Frame(mono)
            self.frames.append(frame)
            fr += 1
            print fr
                
class Frame:
    def __init__(self, image):
        self.mono = image
        
    def show(self):
        plt.imshow(self.mono)
        
if __name__ == '__main__':

    filename = '/home/floris/src/C001H001S0063_compiled_compressed.mov' 
    
    movie = Movie()
    movie.load_using_pyglet(filename, nframes=[0,100])
