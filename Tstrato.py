'''
This code will superficially be very similar to diurnalVariation.py
Its intention is to determine the values of Tstrato with longitude.
It should output scatter plots and linear regression fits for each p[osition selected.
Then, all values of Tstrato will be plotted against longitude
'''
#------------------------------------------------------------------------------
# BOILERPLATE CODE
# The code for setting up the environment, along with other useful variables
#
# importing relevant mathematical packages
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scipy.stats as sc

from pathlib import Path 

# importing the Image and Geo classes to allow the viewing of all the images
from Image import Image
from Geo import Geo
from Position import Position

from spectralWavelength import specWav as sW
from surfaceTemp import getSurfaceTemperature as getST
from surfaceTemp import getTemperatureDeficit as getTD
import Plank as Plank

# makes the plots interactable
plt.ion()
        
# loads in a pre-determined geographical calibration datafile
geo = Geo('geo.txt')

# Defniing local constants for use in the program
channel = 8
histDim = [21,21]
times = list(range(25))

stVars = [9,10,2.9804]

cstr = str(channel)
if channel < 10: cstr = "0"+ cstr
imgDirec = Path("./msg_c" + cstr + ".img")


charWav = sW(channel,Path("./spectralResponses"))

#-----------------------------------------------------------------------------
# SELECTING POSITIONS
# The code to load in a visual image and allow the user to select positions to calibrate our values of Tstrato
#
# loads in the visible image for the middle of the day to pick which locations should be averaged over
img = Image(Path("./msg_c01.img") / "msg_c01_z12.img")
img.disp(1)

# Allows the user to click on the image to select some positions to be averaged at
numPos = int(input("How many positions would you like to average at: "))
positionsST = []
positionsI = []
for i in range(numPos):
    # gets the x and y image positions of a mouse click on the displayed image
    ix,iy = img.click()
    positionsST.append(Position([ix,iy],histDim,geo))
    positionsI.append(Position([ix,iy],histDim,geo))
    print("Longitude: ",positionsST[-1].longitude)
    print("Local time offset: ",positionsST[-1].localTimeOffset)

#-----------------------------------------------------------------------------
# CALCULATING RADIANCES AT EACH TIME STEP
# This code loads in images at each timestamp (hourly) and calculates the Surface radiance B0, as well as records the ozone channel radiance
#

for t in times:
    tstr = str(t)
    if t < 10: tstr = "0" + tstr

    # loads in the ozone channel Image object
    cImg = Image(imgDirec / str("msg_c"+cstr+"_z"+tstr+".img"))
    # also creates the surface temperature image, and converts it into a radiance in the channel spectrum
    ST = getST(stVars[0],stVars[1],stVars[2],t)
    ST.radiance(charWav)

    # for each position, stores the averaged pixel value in each radiance channel
    for i in range(numPos):
        positionsST[i].getVal(ST,t)
        positionsI[i].getVal(cImg,t)

    # alerts user as to progress of program
    print("z"+tstr+ " processed...")

#-----------------------------------------------------------------------------
# CALCULATING REGRESSION FITS
# In this section, we'll calculate the difference in radiance values, and perform regression fits between the data
#
# we will start by defining the array to store the regression fits and scatter data
posData = []

# We'll also create a variable for storing all the difference data in one, for a total linear regression fit
totalData = type('', (), {})()
totalData.dST = []
totalData.dI = []

# for each position, we will caluclate the differences in B0 and I, and the regression fit, and store it in an otherwise empty object
for i in range(numPos):
    # creates a temporary object to store our data: will have the handles .dST .dI and .linReg
    tempOb = type('', (), {})()
    
    dST = []
    dI = []
    # for each data value stored at each position, will caluclate the difference in radiances
    for t in range(len(times)-1):
        dST.append( positionsST[i].vals[t+1] - positionsST[i].vals[t] )
        dI.append( positionsI[i].vals[t+1] - positionsI[i].vals[t] )

    # uses the scipy.stats module to calculate the linear regression between dST and dI
    linReg = sc.linregress(dST,dI)

    # stores data in tempOb, then stores that in posData
    tempOb.dST = dST
    tempOb.dI = dI
    tempOb.linReg = linReg
    tempOb.pos = [positionsST[i].x,positionsST[i].y]

    posData.append(tempOb)

    # appends the data to the totalData object
    totalData.dST += dST 
    totalData.dI += dI

#calculates the linear regression on all the data
totalData.linReg = sc.linregress(totalData.dST,totalData.dI)

#---------------------------------------------------------------------------------
# PLOTTING DATA
# Here we'll make a scatter plot with all the data on one plot, as well as individual plots for each position.
# We'll also plot the value of Tstrato against longitude
#
# Designating window 0 for the scatter plot with all the data included.
fig0 = plt.figure(0)
fig0.clf()
ax0 = fig0.add_subplot(1,1,1)

# for each object in posData, will include data in figure 0, and create new figure to display data
for i in range(numPos):
    data = posData[i]
    # creates figure in new window and plots scatter data on it
    figI = plt.figure(i+1)
    figI.clf()
    #creates axis object to allow for easier positioninng of textbox when plotting
    axI = figI.add_subplot(1, 1, 1)
    axI.scatter(data.dST,data.dI,label=r"$(x,y)=($" + str(data.pos[0]) + "," + str(data.pos[1]) + r"$)$" )

    # creates points to calculate linear regression line, and plots it to same axes
    xmin,xmax = min(data.dST),max(data.dST)
    x = np.linspace(xmin-0.1,xmax+0.1,5)
    y = data.linReg[0]*x + data.linReg[1]
    axI.plot(x,y,label="linear regression")
    plt.xlabel(r"$\Delta B_0$")
    plt.ylabel(r"$\Delta I_8$") # change label for appropriate channel (same with title)
    plt.title(r" Relation of $\Delta I_8 \propto \Delta B_0$ at $(x,y)=($" + str(data.pos[0]) + "," + str(data.pos[1]) + r"$)$")
    plt.legend()
    #### INCLUDE TEXT BOX IN EACH PLOT CONTAINING REGRESSION MODEL DATA
    regString = '\n'.join((
        "Regression Model:",
        r"$\Delta I_8 = \tau_{strato.} \Delta B_0$",
        r'$\tau_{strato.}=%.5f$' % (data.linReg[0], ),
        r'$\mathrm{Intercept}=%.5f$' % (data.linReg[1], ),
        r'$R=%.5f$' % (data.linReg[2], )))
    axI.text(0.65, 0.35, regString, transform=axI.transAxes, fontsize=12,
        verticalalignment='top',bbox=dict(facecolor='orange', alpha=0.5))
    plt.show()

    # return to figure 0 and include scatter
    plt.figure(0)
    plt.scatter(data.dST,data.dI,label=r"$(x,y)=($" + str(data.pos[0]) + "," + str(data.pos[1]) + r"$)$" )

# plots the total linear regression model on figure 0 over the scatter plots
xmin,xmax = min(totalData.dST),max(totalData.dST)
x = np.linspace(xmin-0.1,xmax+0.1,5)
y = data.linReg[0]*x + data.linReg[1]
ax0.plot(x,y,label="linear regression")
regString = '\n'.join((
        "Regression Model:",
        r"$\Delta I_8 = \tau_{strato.} \Delta B_0$",
        r'$\tau_{strato.}=%.5f$' % (data.linReg[0], ),
        r'$\mathrm{Intercept}=%.5f$' % (data.linReg[1], ),
        r'$R=%.5f$' % (data.linReg[2], )))
ax0.text(0.65, 0.35, regString, transform=axI.transAxes, fontsize=12,
         verticalalignment='top',bbox=dict(facecolor='orange', alpha=0.5))

# once all positions accounted for, setup figure titles
plt.xlabel(r"$\Delta B_0$")
plt.ylabel(r"$\Delta I_8$")
plt.title(r"Relation of change in $\Delta I_8 \propto \Delta B_0$")
plt.legend()
plt.show()

# plots Tstrato (slope values) against longitude
plt.figure(numPos+1)
plt.clf()
# for each dataPos object, plots its longitude and slope value
for data in posData:
    lat,long,zen = geo.locate(data.pos[0],data.pos[1])
    print(long)
    plt.scatter(long,data.linReg[0],color='blue')
plt.xlabel(r"Longitude, $\phi$ [$^{\circ}$]")
plt.ylabel(r"Stratospheric Ozone Transmittance, $T_S$")
plt.title("Stratospheric Ozone Transmittance vs. Longitude")
plt.show()

input("Calculate Radiances and Temperatures")
#----------------------------------------------------------------------------------------------
# CALCULATING STRATOSPHERIC RADIANCE AND TEMPERATURE
# This will use the PLank functions to calculate the stratospheric ozone's radiance and temperature above the desert regions
#
# setup for plotting all of them doubled up against local time:
# additional plot for surface radiance on diurnal cycle
figAll,axsAll = plt.subplots(3,sharex=True)
figAll.suptitle(r"$B_{strato.}$ and $T_{strato.}$ vs. local time, ALL")


# for each position, we can go through its I8 and B0 values to calculate the ozone radiance
for i in range(numPos):
    I = positionsI[i].vals
    B0 = positionsST[i].vals
    Bstrato = []
    Tstrato = []
    Ts = posData[i].linReg[0]
    posTimes = (positionsST[i].times)

    # calculates stratospheric brightness for each time (assuming Ts remains constant) and temperature
    for t in times:
        B = (I[t] - B0[t]*Ts)/(1-Ts)
        T = Plank.T(B,charWav)
        Bstrato.append(B)
        Tstrato.append(T)

    #prints out mean and SE on values calculated
    tMean = np.mean(Tstrato)
    tSE = np.std(Tstrato)/np.sqrt(len(times))
    bMean = np.mean(Bstrato)
    bSE = np.std(Bstrato)/np.sqrt(len(times))

    print("T mean: " + str(tMean) + "+-" + str(tSE), "\nB mean: " + str(bMean) + "+-" + str(bSE),"\n")

    # doubling up on each plot (x axis)
    # also doubles up surface radiance so that can be plotted
    plotT = np.linspace(0,49) + positionsST[i].localTimeOffset
    plotTs = [Tstrato[a%25] for a in range(50)]
    plotBs = [Bstrato[a%25] for a in range(50)]
    plotB0 = [Plank.T(B0[a%25],charWav) for a in range(50)]
    
    # creates new figures, each with two subplots (temperature and radiance)
    fig,axs = plt.subplots(2,sharex=True)
    fig.suptitle(r"$B_{strato.}$ and $T_{strato.}$ vs. local time at (x,y)=(" + str(positionsST[i].x) + "," + str(positionsST[i].y) + "), longitude = " + str(geo.locate(positionsST[i].x,positionsST[i].y)[1]))

    # REPLACE plotVars with times,Bstrato,Tstrato
    axs[0].scatter(plotT,plotBs)
    axs[0].set_ylabel(r"$B_{strato}$ [Wm$^{-2}$sr$^{-1}\mu$m$^{-1}$]")
    axs[1].scatter(plotT,plotTs)
    axs[1].set_ylabel(r"$T_{strato}$ [K]")
    axs[1].set_xlabel("local time [Hours]")

    plt.show()

    # also includes scatter on all plot
    axsAll[0].scatter(plotT,plotBs,label="(x,y)=(" + str(positionsST[i].x) + "," + str(positionsST[i].y) + ")")
    axsAll[1].scatter(plotT,plotTs)
    axsAll[2].scatter(plotT,plotB0)




axsAll[0].set_ylabel(r"$B_{strato}$ [Wm$^{-2}$sr$^{-1}\mu$m$^{-1}$]")
axsAll[1].set_ylabel(r"$T_{strato}$ [K]")
axsAll[2].set_ylabel(r"$B_0$ [Wm$^{-2}$sr$^{-1}\mu$m$^{-1}$]")
axsAll[2].set_xlabel("local time [Hours]")
axsAll[0].legend()
plt.show()

plt.show()
    
    
    




input("kill program")



    
