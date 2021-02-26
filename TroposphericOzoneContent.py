'''
This code will calculate the value of x across each image and display it.
It is using the same stratospheric transmittance Ts as the ClearAir script
'''
#------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# The code for setting up the environment, along with other useful variables
#
# importing relevant mathematical packages
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scipy.stats as sc

from pathlib import Path
import copy as copy

# importing the Image and Geo classes to allow the viewing of all the images
from Image import Image
from Geo import Geo

from Position import Position as Position
from Position import loadPositions as loadPositions

from spectralWavelength import specWav as sW
from surfaceTemp import getSurfaceTemperature as getST
from surfaceTemp import getTemperatureDeficit as getTD
import Plank as Plank
from WeightedMean import weightedMean as weightedMean

# makes the plots interactable
plt.ion()
# increases font size on the plots
'''
plt.rc('font',size=18)
plt.rc('axes', titlesize=22)
plt.rc('axes', labelsize=18)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('legend', fontsize=16)
plt.rc('figure', titlesize=22)
'''
        
# loads in a pre-determined geographical calibration datafile
geo = Geo('geo.txt')

# Defniing local constants for use in the program
channel = 8
histDim = [21,21]
times = list(range(25))#[10,11,12,13,14,15,16]#

# reverses list ordering of times for convenience in viewing
times = [times[-(i+1)] for i in range(len(times))]

stVars = [9,10,2.9804]

cstr = str(channel)
if channel < 10: cstr = "0"+ cstr
vImgDirec = Path("./msg_c01.img")
cImgDirec = Path("./msg_c" + cstr + ".img")


charWav = sW(channel,Path("./spectralResponses"))

# The variables obtained form the clean air calibration of the model
Ts = 0.37878 # stratospheric transmittance
dTs = 0.00165 # error on stratospheric transmittance

ts = 249.53916 # stratospheric ozone brightness temperature [Kelvin]
dts = 0.84451 # error on mean of stratospheric ozone temperature [Kelvin]

Bs = Plank.B(ts,charWav) # stratospheric brightness
dBs = Plank.Berr(ts,charWav,dts) # error on stratospheric ozone brightness

dtTrop = -5 # temperature change for tropospheric ozone from surface temperature (taken as positive increase)[Kelvin]


# Loads in the positions from the calibration. As another consitancy check, they should yield x ~ 0
posFile = "pos.txt"
histDim = [21,21]
positions = loadPositions(posFile,geo,histDim)
positionsErr = loadPositions(posFile,geo,histDim)
# Loads in postions for further testing with our model
testPosFile = "testPos.txt"
testPositions = loadPositions(testPosFile,geo,histDim)
testPositionsErr = loadPositions(testPosFile,geo,histDim)


#-------------------------------------------------------------------------------
# CALCULATING TROPOSPHERIC OZONE CONTENT, x
# At each timestamp, we'll open the radiance image and generate the surface temperature image
# We'll then use the formula derived from our model to calculate x.
# We will be using the approximation that tt = t0-5K (as suggested by Anu)
# It will also be displayed along with an image in the visible channel for reference

# before first loop itteration creates variables to store space and zenith mask
spaceMask = []
zenithMask = []

# for each timestamp, we'll open the visible and channel images
for t in times:
    tstr = str(t)
    if t < 10: tstr = "0" + tstr

    # opens the channel radiance image, and gets the surface radiance
    cImg = Image(cImgDirec / str("msg_c" + cstr + "_z" + tstr + ".img"))
    ST = getST(stVars[0],stVars[1],stVars[2],t)

    # from the surface temperature image, creates one with a uniformly colder temperature of dtTrop
    # uses copy module to save having to open image files again
    tropT = copy.deepcopy(ST)
    tropT.data += dtTrop

    # on the first itteration, designates the values for the space and zenith mask
    if t == times[0]:
        # creates a mask that sets the values in the space-region of the image to be zero
        spaceMask = [d > 15 for d in ST.data.flatten()]
        spaceMask = np.array(spaceMask).reshape((ST.ny,ST.nx))

        #creates a mask of all points which the geo object locates as having a zenith angle less than 60deg
        zenithMask = [geo.locate(i%ST.nx,np.floor(i/ST.nx))[2] <= 60 for  i in range(len(ST.data.flatten()))]
        zenithMask = np.array(zenithMask).reshape((ST.nx,ST.ny))

    # sets the surface temperature and reduced temperature to radiance images
    ST.radiance(charWav)
    tropT.radiance(charWav)

    # calculates the relative stratospheric ozone content, x
    # SORT OUT ACTUAL STRATOSPHERIC BRIGHTNESS (placeholder 3.8)

    #stratB = (cImg.data - Ts * ST.data)/(1-Ts)
    x = ( cImg.data - Ts*ST.data - (1-Ts)*Plank.B(ts,charWav) ) / ( (tropT.data - ST.data) *Ts)

    # also calculates the error on x using the derived formula (see writeup)
    # error in I ignored because dI = 0
    errComp = []
    errComp.append(( (1-Ts) * dBs / (tropT.data - ST.data) / Ts )**2) # error in stratospheric radiance
    errComp.append(( x.data * tropT.Berr / (tropT.data - ST.data) )**2) # errror in tropospheric radiance
    errComp.append(( (1-x) * ST.Berr / (tropT.data - ST.data) )**2) # error in surface radiance 
    errComp.append(( (Bs - ST.data) / (tropT.data - ST.data) - x )**2 * ( dTs/Ts )**2) # error in transmittance

    dx = np.sqrt(sum(errComp)) * x.data
    

    # creates the new Image object for displaying x
    xImg = Image(None)
    xImg.title = "Relative Tropospheric Ozone Content at z" + tstr
    xImg.nx = ST.nx
    xImg.ny = ST.ny
    xImg.tem = True # set so can get colour map that is normalised for smaller values
    xImg.data = x*spaceMask + 200*(1-spaceMask)#*spaceMask
    xImg.TOC = True

    # displays xImg
    xImg.disp(t)

    #saves the image to z#.png
    #plt.savefig("z" + str(t) + ".png",bbox_inches="tight",dpi=500)

    #creates new Image object to store error in x. As it won't be viewed, don't have to worry about any other details
    dxImg = Image(None)
    dxImg.nx = ST.nx
    dxImg.ny = ST.ny
    dxImg.data = dx

    # for each position, records the TOC value there
    for pos in positions:
        pos.getVal(xImg,t)
    for pos in positionsErr:
        pos.getVal(dxImg,t)
        
    for pos in testPositions:
        pos.getVal(xImg,t)
    for pos in testPositionsErr:
        pos.getVal(dxImg,t)

    print("z"+ tstr + " processed...")
    #if "q" in input("Next (q to quit): "): break

input("Analyse TOC at specified positions: ")
#------------------------------------------------------------------------------------
# ANALYSIS AT POSITIONS
# this code will go through the positions and visualises TOC at each as function of time

# goes through each of the calibration positions and sees if they have x ~ 0
for i,pos in enumerate(positions):
    x,y = pos.x,pos.y

    plt.figure(i+1,figsize=(19.2,9.83))
    plt.clf()

    plt.subplot(1,1,1)
    # plots values of TOC against time
    plt.errorbar(times,pos.vals,yerr=positionsErr[i].vals,fmt="+")

    #calculates weighted mean of TOC values
    mX, mXErr = weightedMean(pos.vals,positionsErr[i].vals)
    plt.plot([min(times)-0.5,max(times)+0.5],[mX,mX], label=(r'$mean x_{TOC}=%.5f \pm %.5f$' % (mX,mXErr)) )
    plt.fill_between([min(times)-0.5,max(times)+0.5],[mX-mXErr,mX-mXErr],[mX+mXErr,mX+mXErr],alpha=0.3)

    plt.ylabel(r"$x_{TOC}$")
    plt.xlabel("Time, z [hours]")
    plt.legend()
    plt.title("Relative TOC against time at calibration location: (x,y) = (" + str(x) + "," + str(y) + ")")
    plt.show()
    plt.savefig("c" + str(i) + ".png",bbox_inches="tight",dpi=500)

input("Show test positions: ")

# goes through each of the test positions and does the same. Plots these against local time
for i,pos in enumerate(testPositions):
    x,y = pos.x,pos.y

    plt.figure(i+1,figsize=(19.20,9.83))
    plt.clf()

    plt.subplot(1,1,1)
    # plots values of TOC against time
    plt.errorbar(pos.times,pos.vals,yerr=testPositionsErr[i].vals,fmt="+")
    
    plt.ylabel(r"$x_{TOC}$")
    plt.xlabel("local time, t [hours]")

    plt.title("Relative TOC against time at testing location: (x,y) = (" + str(x) + "," + str(y) + ")")
    plt.show()
    plt.savefig("t" + str(i) + ".png",bbox_inches="tight",dpi=500)
    
    
input("kill program")
    
