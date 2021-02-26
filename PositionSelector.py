'''
Code to allow user to select positions to be used for the calibration of the bi-layer atmosphere model
Will allow the user to select as many positions on image as they like.
This can be done by clicking on the image at positions with low-cloud cover.
They will be stored in a file that can then be read in later to maintain consistancy in the model.
'''
#--------------------------------------------------------------
# ENVIRONMENT SETUP
# Importing any relevant modules and setting up any variables for use in the program

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from Image import Image

#the times for the images in the visible channel to be displayed
times = [9,10,12,14,16,17]

# gives Path object to visible image channel
visPath = Path("./msg_c01.img")

# The dimesnsions of the boxes to be drawn on the image to show previously selected regions
hx,hy = 21,21

#make plots interactive
plt.ion()

#--------------------------------------------------------------
# SELECTING POSITIONS
# This code will load in several images in the visible channel, so that diurnal changes in cloud cover can be observed.
# It will then allow the user to select positions in the Z=12 image. Previous selections will be overlayed on the image.

# code to display all images in visible channel 01
for i,t in enumerate(times):
    tstr = str(t)
    if t<10: tstr = "0" + tstr

    img = Image(visPath / str("msg_c01_z" + tstr + ".img"))
    img.disp(i+1)

selectImg = Image(visPath / "msg_c01_z12.img")
selectImg.disp(0)

# user inputs how many positions they want, and selects locations by clicking on window 0
numPos = int(input("How many positions would you like to select? "))
coords = []
boxes = []
for i in range(numPos):
    # gets positions clicked, and appends to coords array
    ix,iy = selectImg.click()
    coords.append([ix,iy])

    # generates box to be drawn on image to demonstrate locations selected
    xmin = ix-np.floor(hx/2)
    ymin = iy-np.floor(hy/2)
    box = {'xmin':xmin, 'ymin':ymin, 'xmax':xmin+hx, 'ymax':ymin+hy, 'color':"red"}
    boxes.append(box)

    selectImg.disp(0,boxes)

#----------------------------------------------------------------
# SAVING FILE
# This will save the positions in the lineformat ix iy\n to a filename of the users choice

filename = input("Please enter full filename (including extension) for positions to be saved to: ")


# generates strings to be written into file
strings = []
for p in coords:
    coordString = str(p[0]) + " " + str(p[1]) + "\n" 
    strings.append(coordString)


f = open(filename,'w')
f.writelines(strings)
f.close()

print("File succesfully closed: " + str(f.closed))
input("Kill program:")
