'''
code for determining the surface temperaturew in a given image
INITIALLY this will use the analysis from original practical, section 5.1, to determine the surface temperature
'''
# importing relevant mathematical packages
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from pathlib import Path

# importing the Image and Geo classes to allow the viewing of all the images
from Image import Image
from Geo import Geo

from spectralWavelength import specWav as sW

#-------------------------------------------------------------------------------------------
def getSurfaceTemperature(i,j,gamma,t):
    '''
    calculates the surface temperature from images in channels i and j based on the formula
    T0 = TI(i) + gamma{i,j}(TI(i) - TI(j))

    also calculates the systematic error on values from calibration uncertainty
    dT0 = dTI(i)_cal + gamma{i,j}( d(TI(i)-TI(j))_WVC )
        ~ 0.9K + gamma * ( ~1K (from WVC graph) )

    INPUTS:
        i [int] : channel i used in caluclation
        j [int] : channel j used in the caulation
        gamma [float] : gamma value based on WVC between different channels
        t [int] : timestamp for image which ST is to be calculated for

    OUTPUTS:
        ST [Image] : new image object whose data should contain surface temperature values
    '''
    # generates empty Image object to return later
    ST = Image(None)

    istr = str(i)
    jstr = str(j)
    tstr = str(t)
    if i < 10: istr = "0"+istr
    if j < 10: jstr = "0"+jstr
    if t < 10: tstr = "0"+tstr

    # calculates characteristic wavelengths for conversion into brightness temperatures
    wavI = sW(i,Path("./spectralResponses"))
    wavJ = sW(j,Path("./spectralResponses"))
    
    # opens the i and j channel images
    imgI = Image(Path("./msg_c" + istr + ".img") / str("msg_c"+istr+"_z" + tstr +".img"))
    imgJ = Image(Path("./msg_c" + jstr + ".img") / str("msg_c"+jstr+"_z" + tstr +".img"))
    imgI.bright(wavI)
    imgJ.bright(wavJ)

    # uses formula to calculate new data from existing images' data
    newData = imgI.data + gamma*(imgI.data - imgJ.data)

    # calculates sytematic uncertainty on ST data
    dT = 0.9 + gamma*(1)
    
    # set the Image objects variables
    ST.tem = True
    ST.nx = imgI.nx
    ST.ny = imgI.ny
    ST.data = newData
    ST.title = "Surface Temperature ("+ istr + "," + jstr+")"
    ST.systematicErrorTemp = dT
    return ST

#---------------------------------------------------------------------------------------------
def getTemperatureDeficit(i,ST,t):
    '''
    A function for determining the temperature deficit in a channel given a surface temperature plot

    INPUTS:
        i [int] : channel for which temperature deficit is to be calculated
        ST [Image] : Image object for the surface temperature
        t [int] : timestamp for the time of the image to be used

    OUTPUTS:
        TD [Image] : Image object with the temperature deficits calculated
    '''
    # generates the empty image object to be returned
    TD = Image(None)

    istr = str(i)
    tstr = str(t)
    if i < 10: istr = "0" + istr
    if t < 10: tstr = "0" + tstr
    # calculates the characteristic wavelength for conversion into temperature defecit
    wavI = sW(i,Path("./spectralResponses"))

    # opens the image in the i channel and converts to a brightness temperature
    imgI = Image(Path("./msg_c" + istr + ".img") / str("msg_c" + istr + "_z" + tstr + ".img"))
    imgI.bright(wavI)

    # calculates positive temperature deficit for channel i
    newData = ST.data - imgI.data

    # sets the ST object's parameters
    TD.tem = True
    TD.temDef = True
    TD.nx = imgI.nx
    TD.ny = imgI.ny
    TD.data = newData
    TD.title = "Temperature deficit in c"+istr
    return TD

#----------------------------------------------------------------------------------------------

'''
Test code:
st = getSurfaceTemperature(9,10,2.9804,12)
st.disp(1)
imgI = Image(Path("./msg_c08.img") / str("msg_c08_z12.img"))
imgI.bright(sW(8,Path("./spectralResponses")))
imgI.disp(2)
td = getTemperatureDeficit(8,st,12)
td.disp(3)
input("Kill Program")
'''

    

