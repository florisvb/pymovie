import pickle
from pygmovie import npMovie
from pygmovie import uFrame
from pygmovie import Object
from pygmovie import Frame

filename = '/home/floris/data/windtunnel/SA1/sa1_npmovie_20101101_C001H001S0019'

fname = (filename)  
fd = open( fname, mode='r' )
npmovie = pickle.load(fd)
fd.close()

