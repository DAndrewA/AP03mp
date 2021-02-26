'''
This code will allow the user to calculate values of Ts for pre-selected positions.
For each position, we will do a linear regression on dST ~ dI8.
We will then calculate a weighted mean of the Ts values to get a globally applicable Ts value
This can then be used to caluclate the stratospheric ozone brightness temperature ts
We can compare this in our calibration regions to see if they provide temperatures that agree 
(Literature estimate, ts ~ 250K)
'''

# defines the equation of a line: used when performing the linear regression fitting
def line(x,a,b):
    return a*x + b

#------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# The code for setting up the environment, along with other useful variables
# Also has the code for changing matplotlib parameters (i.e fontsize, etc)
#
# importing relevant mathematical packages
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scipy.stats as sc
from scipy.optimize import curve_fit as sccf

from pathlib import Path 

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
plt.rc('font',size=18)
plt.rc('axes', titlesize=22)
plt.rc('axes', labelsize=18)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('legend', fontsize=16)
plt.rc('figure', titlesize=22) 

        
# loads in a pre-determined geographical calibration datafile
geo = Geo('geo.txt')

posFile = "pos.txt"

# Defniing local constants for use in the program
channel = 8
histDim = [21,21]
times = list(range(25))

# variables used in calculating surface temperature
stVars = [9,10,2.9804]

cstr = str(channel)
if channel < 10: cstr = "0"+ cstr
cImgDirec = Path("./msg_c" + cstr + ".img")

# the characteristic wavelength for the channel
charWav = sW(channel,Path("./spectralResponses"))

# The variables used in the clean-air two-layer atmosphere model
#Ts = 0.36

#-----------------------------------------------------------------------------
# LOADING IN POSITIONS
# loads in the positions from the given file (pos.txt) into two seperate arrays, for each image dataset
positionsST = loadPositions(posFile,geo,histDim)
positionsI = loadPositions(posFile,geo,histDim)
# also loads in an additional array for calculating errors in surface radiance
positionsSTErr = loadPositions(posFile,geo,histDim)

#-----------------------------------------------------------------------------
# LOADING RADIANCE DATA
# for each timestamp, will get the surface radiance and ozone channel radiance and store their values at the given Positions

for t in times:
    tstr = str(t)
    if t<10: tstr = "0" + tstr

    #loads in the ozone radiance image
    ozImg = Image(cImgDirec / str("msg_c"+cstr+"_z"+tstr+".img"))

    # creates the surface temperature image, and converts it into a surface radiance
    ST = getST(stVars[0],stVars[1],stVars[2],t)
    ST.radiance(charWav)

    # creates new image object to store brightness error values in main image data
    errImg = Image(None)
    errImg.nx = ST.nx
    errImg.ny = ST.ny
    errImg.data = ST.Berr
    ST.title = "Surface Radiance Error (9,10)"
    


    # for each position in positionsST, gets value of surface radiance at that position
    for pos in positionsST:
        pos.getVal(ST,t)
        
    # for each position in positionsI, gets value of ozone radiance at that position
    for pos in positionsI:
        pos.getVal(ozImg,t)

    # for each position in positionsSTErr, gets value of brightness error at that position
    for pos in positionsSTErr:
        pos.getVal(errImg,t)

    print("z"+ tstr + " processed...")

#-----------------------------------------------------------------------------
# CALCULATING REGRESSION MODELS
# Will use the scipy curve fitting with errors on values to determine values for Ts and intercept
# Can also use scipy stats model to calculate regression model on the data itself.

# We'll create an array to store objects holding the data we want to consider at each position
posData = []

# for each position, will go through and calculate differences in ST (actually B0) and I at each timestamp
# will also calculate error at each step in B0
for i,p in enumerate(positionsST):
    # creates a temporary object to store our data: we can define any handles we like for it
    tempOb = type('', (), {})()

    dSTvals = []
    dIvals = []
    dSTErrVals = []
    # for each data value, will caluclate the radiance change at that poisition (units radiance/hour)
    for t in range(len(times)-1):
        dST = p.vals[t+1] - p.vals[t]
        dI = positionsI[i].vals[t+1] - positionsI[i].vals[t]
        dSTvals.append(dST)
        dIvals.append(dI)

        # also calculates the error in each data point for ST
        # systematic error + quadrature(error on subtracted points)
        dSTErr = positionsSTErr[i].vals[t] + dST * ( (p.valsError[t+1]/p.vals[t+1])**2 + (p.valsError[t]/p.vals[t])**2 )**0.5
        dSTErrVals.append(dSTErr)

    # error in surface temperature is the larger of the error values here so we'll consider it to be the primary source of error. 

    # ASSUMING HAVE ERROR IN ST ALREADY as errST
    
    # calculates the line fit with the errors for dST = (1/Ts) dI
    pOpt, pCov = sccf(line,dIvals,dSTvals,sigma=dSTErrVals)

    # calculates the Ts value with error
    Tstrato = 1/pOpt[0]
    TstratoErr = Tstrato**2 * pCov[0,0]**0.5

    #uses scipy module to caluclate linear regression on dI ~ dST (other way around)
    linReg = sc.linregress(dSTvals,dIvals)

    # stores data in tempOb
    tempOb.linReg = linReg # linear regression model
    tempOb.dST = dSTvals # dST values for scatter plot
    tempOb.dSTErr = dSTErrVals # error values on surface brightness
    tempOb.dI = dIvals # dI values for scatter plot
    tempOb.Tstrato = [Tstrato,TstratoErr] # For comparing to linReg values. If agree, will use these ones
    tempOb.pos = [p.x,p.y] # image coordinates of position, for ease of labelling

    # appends this object to our posData array
    posData.append(tempOb)

#-----------------------------------------------------------------------------
# PLOTTING DATA / VISUALISING RESULTS
# This section will be responsible for generating any plots required to visualise the data, as well as outputting values so I can record them

# designates window 0 for scatter plot of all values.
fig0 = plt.figure(0,figsize=(19.20,9.83))
fig0.clf()

# for each object in posData (1 to 1 with Positions) will create scatter plot with regression line on
for i,pos in enumerate(posData):
    # creates figure in unique window, and axis for plotting on
    figI = plt.figure(i+1,figsize=(19.20,9.83))
    figI.clf()
    axI = figI.add_subplot(1,1,1)

    # extract data from object pos to be used
    dST = pos.dST
    err = pos.dSTErr
    dI = pos.dI
    x,y = pos.pos[0],pos.pos[1]
    
    # plots scatter plot with the dST error on axis
    axI.errorbar(dST, dI, xerr=err, fmt="+",label=r"$(x,y)=($" + str(x) + "," + str(y) + r"$)$")
    # creates points to calculate linear regression line, and plots it to same axes
    xmin,xmax = min(pos.dST),max(pos.dST)
    lx = np.linspace(xmin-0.1,xmax+0.1,5)
    ly = pos.linReg[0]*lx + pos.linReg[1]
    axI.plot(lx,ly,label="linear regression")
    # labelling of axis
    plt.xlabel(r"$\Delta B_0 [Wm^{-2}sr^{-1}\mu m^{-1} (hour)^{-1}]$")
    plt.ylabel(r"$\Delta I_8 [Wm^{-2}sr^{-1}\mu m^{-1} (hour)^{-1}]$") # change label for appropriate channel (same with title)
    plt.title(r" Relation of $\Delta I_8 \propto \Delta B_0$ at $(x,y)=($" + str(x) + "," + str(y) + r"$)$")
    plt.legend()

    #### INCLUDE TEXT BOX IN EACH PLOT CONTAINING REGRESSION MODEL DATA
    regString = '\n'.join((
        "Regression Model:",
        r"$\Delta I_8 = \tau_{strato.} \Delta B_0$",
        r'$\tau_{strato.}=%.5f$' % (pos.linReg[0], ),
        r"$scipy \delta\tau_{strato.} = %.5f$" % (pos.linReg[-1]),
        r'$\mathrm{Intercept}=%.5f$' % (pos.linReg[1], ),
        r'$R=%.5f$' % (pos.linReg[2], ),
        r'$Propogated \tau_{strato.}=%.5f \pm %.5f$' % (pos.Tstrato[0],pos.Tstrato[1]) ))
    axI.text(0.65, 0.35, regString, transform=axI.transAxes,
        verticalalignment='top',bbox=dict(facecolor='orange', alpha=0.5))

    plt.show()
    
    # saves the plot to s#.png. UNCOMMENT TO SAVE FIGURES AUTOMATICALLY
    #plt.savefig("s" + str(i+1) +".png",bbox_inches='tight',dpi=500)

    # include in scatter in window 0
    plt.figure(0)
    plt.scatter(dST,dI,label=r"$(x,y)=($" + str(x) + "," + str(y) + r"$)$")


# once all positions accounted for, setup figure 0 labels
plt.xlabel(r"$\Delta B_0 [Wm^{-2}sr^{-1}\mu m^{-1} (hour)^{-1}]$")
plt.ylabel(r"$\Delta I_8 [Wm^{-2}sr^{-1}\mu m^{-1} (hour)^{-1}]$")
plt.title(r"Relation of change in $\Delta I_8 \propto \Delta B_0$")
plt.legend()
plt.show()
#plt.savefig("totalScatter.png",bbox_inches='tight',dpi=500)

input("Plot Tstrato values: ")

# plots all the calculated Tstrato values on a figure
plt.clf()
for i,pos in enumerate(posData):
    Tstrato = pos.Tstrato[0]
    err = pos.Tstrato[1]
    x,y = pos.pos[0],pos.pos[1]
    plt.errorbar(i,Tstrato,yerr=err,fmt="+")


# uses the weighted mean function to calculate an average Tstrato value with error on
TstratoArr = [p.Tstrato[0] for p in posData]
TstratoErrArr = [p.Tstrato[1] for p in posData]
mTstrato,mTstratoErr = weightedMean(TstratoArr,TstratoErrArr)

plt.errorbar(len(posData),mTstrato,yerr=mTstratoErr,fmt="+",label=(r'$mean \tau_{strato.}=%.5f \pm %.5f$' % (mTstrato,mTstratoErr)))
plt.fill_between([-1,len(TstratoArr)+1],[mTstrato-mTstratoErr,mTstrato-mTstratoErr],[mTstrato+mTstratoErr,mTstrato+mTstratoErr],alpha=0.3)

plt.title(r"$T_{strato.}$ values calculated at each position")
plt.ylabel(r"$T_{strato.}$")
plt.legend()
plt.show()
#plt.savefig("TstratoBoxPlot.png",bbox_inches='tight',dpi=500)

input("Calculate ozone brightnesses and temperatures: ")
#----------------------------------------------------------------------------
# CALCULATING TROPOSPHERIC OZONE TEMPERATURE
# This section will use the mean Tstrato value to calculate the stratospheric ozone temperature.
# We can then compare it for each positions to literature values, and see if they all agree.

# will designate window 0 for plotting brightnesses, -1 for temperatures
plt.figure(-1,figsize=(19.2,9.83))
plt.clf()
plt.figure(0,figsize=(19.2,9.83))
plt.clf()

# creates arrays to store mean temperature and error on each position. These can be used to determine outliers
meanT = []
meanTErr = []

# for each Position selected, we have the surface brightness (positionsST) and ozone radiance (positionsI) recorded
# for each positions, for each timestamp we will calculate the tropospheric ozone brightness
for i,posST in enumerate(positionsST):
    # These will hold the hourly brightness values and their errors
    stratB = []
    stratBErr = []

    stratT = []
    stratTErr = []

    # for each timestamp calculates stratospheric brightness, B, under clean air assumptions
    for t in range(len(times)):
        B = ( positionsI[i].vals[t] - mTstrato * posST.vals[t] ) / ( 1 - mTstrato )
        stratB.append(B)
        # calculates the error in Bstrato. (ASSUMING ERROR IN I IS 0)
        BErr = np.sqrt( (mTstrato * positionsSTErr[i].vals[t] )**2 + ( ( B + posST.vals[t] )* mTstratoErr )**2 ) / ( 1 - mTstrato)
        stratBErr.append(BErr)

        # Having calculated brightness, can convert into temperatures using functions in Plank module
        T = Plank.T(B,charWav)
        TErr = Plank.Terr(B,charWav,BErr)
        stratT.append(T)
        stratTErr.append(TErr)

    # calculates mean temperature and error for position, and appends them to total arrays
    mT,mTErr = weightedMean(stratT,stratTErr)
    meanT.append(mT)
    meanTErr.append(mTErr)

    x,y = posST.x,posST.y

    # plots variations on own figure in window(i) Now two subplots, for brightness and temperature
    plt.figure(i+1,figsize=(19.2,15))
    plt.clf()

    plt.subplot(2,1,1)
    plt.errorbar(times,stratB,yerr=stratBErr,fmt="+")
    plt.ylabel(r"$B_{strato.} [Wm^{-2}sr^{-1}\mu m^{-1}]$")
    plt.title("Stratospheric ozone brightness and brightness temperature against time at (x,y) = " + str(x) + "," + str(y) + ")")

    plt.subplot(2,1,2)
    plt.errorbar(times,stratT,yerr=stratTErr,fmt="+")
    # on lower subplot, plots average temperature and error on it as fill
    plt.plot([-0.5,len(times)-0.5],[mT,mT],label=(r'$mean T_{strato.}=%.5f \pm %.5f$' % (mT,mTErr)))
    plt.fill_between([-0.5,len(times)-0.5],[mT-mTErr,mT-mTErr],[mT+mTErr,mT+mTErr],alpha=0.3)
    plt.legend()
    plt.ylabel(r"$T_{strato.} [K]$")
    plt.xlabel("Time [hours]")

    # saves plot to tb#.png
    #plt.savefig("tb" + str(i+1) + ".png",bbox_inches='tight',dpi=500)
    
    
    plt.figure(0)
    # plots the brightness values as a function of local time for each position on figure 0
    plt.errorbar(posST.times,stratB,yerr=stratBErr,fmt="+",label=r"$(x,y)=($" + str(x) + "," + str(y) + r"$)$")
    # on figure 1, plots the brightness temperatures as a function of local time
    plt.figure(-1)
    plt.errorbar(posST.times,stratT,yerr=stratTErr,fmt="+",label=r"$(x,y)=($" + str(x) + "," + str(y) + r"$)$")
    plt.figure(0)

# calculates the weighted mean of the mean temperatures, and its error. Plots with temperatures further than a given sigma will be considered outliers
mT,mTErr = weightedMean(meanT,meanTErr)
mTErr = max(meanTErr)

# formatting figure 0
plt.legend()
plt.xlabel("Local time [hours]")
plt.ylabel(r"$B_{strato.} [Wm^{-2}sr^{-1}\mu m^{-1}]$")
plt.title("Calculated stratospheric ozone brightness against local time")

#plt.savefig("bTotal",bbox_inches='tight',dpi=500)

# formatting figure 1
plt.figure(-1)

plt.plot([-0.5,len(times)-0.5],[mT,mT],label=(r'$mean T_{strato.}=%.5f \pm %.5f$' % (mT,mTErr)))
plt.fill_between([-0.5,len(times)-0.5],[mT-mTErr,mT-mTErr],[mT+mTErr,mT+mTErr],alpha=0.3)

plt.legend()
plt.xlabel("Local time [hours]")
plt.ylabel(r"$T_{strato.} [K]$")
plt.title("Stratospheric ozone brightness temperature against local time")

#plt.savefig("tTotal",bbox_inches='tight',dpi=500)

plt.show()

input("Display mean temperature plot: ")

#------------------------------------------------------------------------------
# SELECTING OUTLIERS
# This section will allow the user to select outliers from the calibration process based on how far their temperature falls from the calculataed mean
# It will also display all the required information about the positions in text form for use in the writeup

plt.figure(0)
plt.clf()

for i in range(len(meanT)):
    plt.errorbar(i+1,meanT[i],yerr=meanTErr[i],fmt="b+")

plt.plot([0,len(meanT)+1],[mT,mT],label=(r'$mean T_{strato.}=%.5f \pm %.5f$' % (mT,mTErr)))
plt.fill_between([0,len(meanT)+1],[mT-mTErr,mT-mTErr],[mT+mTErr,mT+mTErr],alpha=0.3)

plt.ylabel(r"$T_{strato.} [K]$")
plt.legend()
plt.title("Stratospheric brightness temperatures for each position")
plt.show()

#plt.savefig("meanTemp",bbox_inches='tight',dpi=500)

# for each position, will print out:
print("#----------------------------")
print("# Positions Information \n\n")
print("x, y, Tstrato, error Tstrato, stratospheric temperature [K], error [K], sigma from mean \n")

for i, posST in enumerate(positionsST):
    dSigma = (meanT[i] - mT) / mTErr
    print(posST.x,posST.y,posData[i].Tstrato[0],posData[i].Tstrato[1],meanT[i],meanTErr[i],dSigma)


input("kill program: ")
