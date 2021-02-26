'''
Taken from the example code, this class handles the displaying of images

I have also included click method from Image_class from main practical
'''

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


from pathlib import Path
import Plank as Plank

#-------------------------------------------------------------------------------
class Image:  
  """ Image data and methods 

  DATA
    tem   boo : True=brightness temperature image, False=radiance image
    nx    int : No of horizontal pixels in images
    ny    int : No of vertical pixels in images
    ixoff int : pixel offset of left edge of box
    iyoff int : pixel offset of bottom edge of box
    data  flt : np.(ny,nx) image data as 2D array
    title str : image title 

  METHODS 
    __init__ : Initialise new Img object
        disp : Display full disk image
        clip : Create new image from subset of original image
      bright : Convert radiance image to Brightness Temperature image

  HISTORY
      v16Oct20 : AD Original version
  """

  def __init__(self,imgfil,title=None):
    """ Initialise new image object 

    PARAMETERS
      imgfil str : file name for image
      title  str : description of image

    DESCRIPTION
      Read in image data as 2D array of floating values.
      If called with imgfil==None this just creates a new 'empty' object 
      If called without title parameter, the title will be taken from the file
    """

    self.tem = False    # Assume any new image is radiance
    self.temDef = False # Not a temperature deficit until it is
    self.TOC = False    # Not a tropospheric ozone content image unless it is
    self.ixoff = 0
    self.iyoff = 0
    if imgfil is None:  # Create an empty object
      self.nx = 0
      self.ny = 0
      self.data = 0
      self.title = ''
    else:               # Read in image data from file
      with open(imgfil) as f:
        imgtitle = f.readline()
        if title is None: self.title=imgtitle
        else: self.title=title
        nx, ny = np.fromfile(f,dtype=int,count=2,sep=" ")  
        self.nx = nx
        self.ny = ny
        imgdata = np.fromfile(f, dtype=float, count=nx*ny, sep=" ")
        self.data = imgdata.reshape((ny,nx))

  def disp(self,window=1,box=None):
    """ Display the image 

    PARAMETERS
      window int : plot window for display
      box        : box coordinates, dictionary containing keys
        'xmin'  int : Pixel# of left edge of box
        'xmax'  int : Pixel# of right edge of box
        'ymin'  int : Pixel# of bottom edge of box
        'ymax'  int : Pixel# of top edge of box
        'color' str : (optional, default='white') color of box

    DESCRIPTION
      Basic function is to display the stored image self.data
      Can also superimpose a rectangular box defined by keys 'xm

    """
    plt.figure(window)
    plt.clf()
    plt.axis('off')                  # no axes required
    plt.title(self.title)
    if self.tem:                     # colour scale for bright.tem image
      colours = [(1,1,1),(1,0,1),(0,0,1),(0,1,1),(0,1,0),(1,1,0),
                 (1,0,0),(0,0,0)]
      # changes normalisation based on whether temperature deficit or not
      if self.temDef: tnorm = plt.Normalize(-20,50)
      else: tnorm = plt.Normalize(230,320) # range 230-320K
      # Only for use in Tropospheric Ozone Content images
      if self.TOC: tnorm = plt.Normalize(-5,25)
      
      tcmap = LinearSegmentedColormap.from_list('tem_colours', colours)
      plt.imshow(self.data, origin='lower', cmap=tcmap, norm=tnorm)
      plt.colorbar(label='Temperature [K]')
    else:                            # use gray scale for radiance image
      plt.imshow(self.data, origin='lower', cmap='gray')
      plt.colorbar(label='Radiance')
    if box is not None:
      for b in box:
        xmin = b['xmin'] - self.ixoff
        xmax = b['xmax'] - self.ixoff
        ymin = b['ymin'] - self.iyoff
        ymax = b['ymax'] - self.iyoff
        if 'color' in b: col=b['color']
        else: col='white'
        plt.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color=col)
           
    plt.tight_layout(pad=0.05)       # maximise size of image within window
    plt.show()

  def clip ( self, box, title=None ):
    """ Create a new image from a subset of another image

    PARAMETERS
      box   dict : Subsect coordinates, as defined in self.disp
      title str  : Title for new image (else copy original image title)

    RETURNS
      newimg : Image object

    DESCRIPTION
      Creates a new image object from a rectangular subset of an existing image
      Note that pixel numbers from the original image have to be preserved, via
      self.ixoff and self.iyoff, in order for the geolocation to work on the
      subset image
    """
    xmin = box['xmin']
    xmax = box['xmax']
    ymin = box['ymin']
    ymax = box['ymax']
    newimg = Image(None)  
    newimg.nx = xmax-xmin 
    newimg.ny = ymax-ymin
    newimg.ixoff = xmin
    newimg.iyoff = ymin
    y1 = ymin - self.iyoff
    y2 = ymax - self.iyoff
    x1 = xmin - self.ixoff
    x2 = xmax - self.ixoff
    newimg.data = self.data[y1:y2,x1:x2]
    if title is None: newimg.title = self.title
    else: newimg.title = title
    return newimg

  def bright ( self, wavelength ):
    """ Convert radiance image to Brightness Temperature image

    PARAMETERS
      wavelength flt : characteristic wavelength [microns] for channel

    DESCRIPTION
      Convert image from radiance [W/(m2 sr um) to brightness temperature [K]
      using the inverse of the Planck function at the characteristic 
      wavelength. 
    """
    # Local constants
    H = 6.63e-34       # Planck constant       [m2 kg / s]
    C = 3.00e8         # Speed of light        [m / s]
    K = 1.38e-23       # Boltzmann constant    [m2 kg /s2 /K]
    R1 = H * C / K     # Intermediate Constant [m K]
    R2 = 2 * H * C**2  # Intermediate Constant [m4 kg / s3]
    if self.tem:
      print('Apparently already a brightness temperature image')
    else: 
      w = wavelength * 1.0e-6  # convert microns to metres
      self.data = R1 / w / np.log( 1.0 + R2/(w**5 * self.data*1e6) )
      self.tem = True

  def radiance(self,wavelength):
    '''
    Converts from a brightness temperature to a radiance value
    also calculates error in brightness as each pixel (if applicable)
    
    INPUTS:
      wavelength [float] : characteristic wavelength for the channel in microns
    '''
    # Local constants
    H = 6.63e-34       # Planck constant       [m2 kg / s]
    C = 3.00e8         # Speed of light        [m / s]
    K = 1.38e-23       # Boltzmann constant    [m2 kg /s2 /K]
    R1 = H * C / K     # Intermediate Constant [m K]
    R2 = 2 * H * C**2  # Intermediate Constant [m4 kg / s3]

    if not self.tem:
      print("Already a radiance image")
    else:
      w = wavelength * 1.0e-6 # converts the wavelength into microns
      b = R2 / w**5 /1e6 / (np.exp(R1 / (self.data * w))-1)
      self.tem = False
      # if the image was created through getST, will use error on ST to caluclate error in B
      if hasattr(self,"systematicErrorTemp"):
        err = b*R1 / w / (self.data)**2 / (1-np.exp( -R1 / (w * self.data) )) * self.systematicErrorTemp
        self.Berr = err

      self.data = b


  def click(self):
    """ User click on image to move cursor box

    DESCRIPTION
      Uses pyplot ginput module assuming a click on mouse left button. See the
      on-line documentation for other button options. Also requires display to
      be set to the full disk image.
    """
    print("left click on image...")
    self.disp('')
    pts = plt.ginput()[0]    # select first pair from list
    ix = int(round(pts[0]))
    iy = int(round(pts[1]))
    return (ix,iy)
