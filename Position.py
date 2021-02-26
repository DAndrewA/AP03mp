# importing relevant mathematical packages
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from pathlib import Path 

class Position:
    '''
    Stores a position on the image, and allows for the calculation of its diurnal variation

    DATA
        position = [x,y] [int array] : position cordinates of the pixel in the image
        longitude [float] : longitude of the pixel selected
        localTimeOffset [float] : offset of local time as measured from 0degE
        hx,hy [int] : dimensions of histogram box applied around pixel

        vals [float array] : values of pixel taken at location
        valsError [float array] : values of standard error calculated for mean pixel values
        times [int array] : the local time the pixel value is recorded at

        box [dictionary] : the box dictionary used to define the histogram area

    METHODS
        __init__ : Initialise position object
        getVal : given an image, caluclates an average value for that pixel
        
    '''
    def __init__(self,pos,hist,geo,colour = None):
        '''
        Initialises the object.

        Inputs:
            pos = [x,y] : image coordinates for pixel
            hist = [hx,hy] : dimensions for histogram which pixel value will be averaged for
            geo [Geo] : the Geo object used to determine the longitude
            colour [string](optional) : the colour for the box used if it is to be drawn
        '''
        self.x = pos[0]
        self.y = pos[1]

        try:
            lat,lon,zen = geo.locate(pos[0],pos[1])
            self.longitude = lon
        except:
            self.longitude = 0
        self.localTimeOffset = -lon/90 * 6 # calculates the local time offset in hours

        self.hx = hist[0]
        self.hy = hist[1]

        self.vals = []
        self.valsError = []
        self.times = []

        xmin = int(self.x - np.floor(self.hx/2))
        ymin = int(self.y - np.floor(self.hy/2))
        self.box = {'xmin':xmin, 'ymin':ymin, 'xmax':xmin+self.hx, 'ymax':ymin+self.hy}
        if colour is not None: slef.box['color'] = colour
        

    def getVal(self,img,t):
        '''
        Method to obtain the pixel value at the location in the image.

        Inputs:
            img [Image] : the image class thats being looked at
            t [int] : the time the image was taken (as caluclated at 0degE)

        Outputs:
            mean [float]: mean pixel value
            SE [float] : standard error on mean
        '''
        # selects the histogram area in the image and returns as new Image object
        boxImg = img.clip(self.box)

        # flattens the image data into a 1d array for ease of manipulation, then calculates mean and standard deviation
        histData = boxImg.data.flatten()
        mean = np.mean(histData)
        sd = np.std(histData)

        # calculates standard error on mean as sd/sqrt(N)
        sE = sd/np.sqrt(self.hx*self.hy)
        # appends calculated values to lists
        self.vals.append(mean)
        self.valsError.append(sE)
        self.times.append((t + self.localTimeOffset))
        return mean,sE


#-------------------------------------------------------------------------------------------------------------------------
# LOAD FUNCTION
# This function will be for loading positions from a formatted text file into an array.
#

def loadPositions(filename,geo,hist):
    '''
    Loads in positions from a text file and returns an array of Position objects

    INPUTS:
        filename [string]: filename for stored positions, with .txt extension
        geo [Geo]: the Geo object used to calibrate each position to the images
        hist [1x2 int array]: dimensions of histogram regions used by positions in collecting data

    OUTPUTS:
        positions [Position array]: array of Position objects specified from file
    '''
    # in try statement incase file doesn't exist
    positions = []
    try:
        f = open(filename,'r')
        data = f.read()
        f.close()
        # splits file contents by line, from which ix and iy can be extracted
        coordStrings = data.split("\n")

        for s in coordStrings:
            s = s.split(" ")
            if s != ['']:
                coords = [int(s[0]),int(s[1])]
                pos = Position(coords, hist, geo)
                positions.append(pos)
        
    except Exception as inst:
        print("Loading positions failed.")
        print(type(inst))    # the exception instance
        print(inst.args)     # arguments stored in .args
        print(inst)
        return None

    return positions



