'''
The class taken from the example code that handles the conversion between relevant coordinate systems
'''

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ------------------------------------------------------------------------------
class Geo:
  """ Geometric calibration data and methods

  DATA 
    cal    boo : True = Geoetric calibration set
    alpha  flt : y/elevation scale factor
    beta   flt : x/azimuth   scale factor
    x0     flt : x-coordinate of sub-satellite point
    y0     flt : y-coordinate of sub-satellite point
    geofil str : Name of file containing GeoCal data (eg 'geo.txt')

 METHODS
    __init__ : Initialise new Geo object   
      menu   : Screen menu for GeoCal part of practical
      save   : Write GeoCal data to file     
      input  : Read user-input of GeoCal parameters 
      angles : Print Elevation,Azimuth angles for givem Lat,Lon
      locang : Convert ele,azi angles to lat,lon,zen angles
      locate : Convert ix,iy coords to lat,lon,zen angles
      satang : Convert lat,lon angles to ele,azi,zen angles

  USAGE
    Called once at the start to initialise a Geo object

  HISTORY
    v16Oct20 : AD Original version
  """
 
  def __init__(self,geofil):
    """ Initialise new Geo object

    PARAMETERS
      geofil str : name of file containing geo.cal data, eg 'geo.txt' 

    DESCRIPTION
      If the Geo Cal data file exists, the data are read on initialisation
      and geo.cal set True
    """

    # Local constants
    self.DIST = 42260.0      # Radial dist [km] of sat. from centre of earth
    self.REARTH  = 6371.0    # Earth radius [km]

    try:                     # if file already exists ...
      f = open(geofil,"r")
      rec = f.readline()  
      rec = f.readline()
      flds = rec.split()
      self.y0    = float(flds[0])   # y-coordinate of sub-satellite point
      self.x0    = float(flds[1])   # x-coordinate of sub-satellite point
      self.alpha = float(flds[2])   # y/elevation scale factor
      self.beta  = float(flds[3])   # x/azimuth   scale factor
      f.close()
      self.cal = True               # Flag for GeoCal data set
      print(" *** GeoCal data loaded from file: " + geofil)
    except:                         # file doesn't exist or can't be read
      print(" *** GeoCal data file not found/read: " + geofil)

  def locang(self,ele,azi):
    """ Convert ele,azi angles to lat,lon,zen angles

    PARAMETERS
      ele flt : Elevation angle [deg]
      azi flt : Azimuth angle [deg]

    RETURNS
      lat flt : Latitude [deg N]
      lon flt : Longitude [deg E]
      zen flt : Zenith angle [deg]

    DESCRIPTION
      The inverse of SATANG.
      Uses spherical coordinate geometry to find the point of intersection of a 
      ray leaving the satellite at particular ele,azi angle with the earth 
      surface
      If no intersection, returns (np.nan,np.nan,np.nan) instead.
    """
    rele     = math.radians(ele)
    sinele   = math.sin(rele)
    cosele   = math.cos(rele)
    razi     = math.radians(azi)
    sinazi   = math.sin(razi)
    cosazi   = math.cos(razi)
    # Distance of plane of intersection from centre of earth
    h = self.DIST * sinele 
    if abs(h) > self.REARTH: return (np.nan,np.nan,np.nan)  # no sfc intersect.
    r1 = math.sqrt( self.REARTH**2 - h**2 ) # Radius of circle of intersection
    d1 = self.DIST * cosele
    if abs ( d1 * sinazi ) > r1: return (np.nan,np.nan,np.nan) # No intersection
    # Distance of line of sight
    x = d1 * cosazi - math.sqrt( r1**2 - d1**2 * sinazi**2 )
    # Distance from pixel to point of intersection of earth's vertical axis with
    # plane of intersection
    d2 = self.DIST / cosele
    y = x**2 + d2**2 - 2 * x * d2 * cosazi
    if y < 0.0: y = 0.0
    y = math.sqrt(y)
    h1 = self.DIST * math.tan(rele)
    if abs(h1) > 1.0e-10:     
      gamma = math.acos( ( self.REARTH**2 + h1**2 - y**2 ) / 
                         ( 2.0 * self.REARTH * h1 ) ) 
    else:
      gamma = math.pi / 2.0 - h1 / ( 2.0 * self.REARTH )
    rlat = math.pi / 2.0 - gamma
    gamma1 = math.asin ( sinazi * x / y )
    rlon = math.atan ( math.sin(gamma1) / ( math.cos(gamma1) * cosele ) )
    rzen = math.acos ( cosazi * cosele ) + \
           math.acos ( math.cos(rlat) * math.cos(rlon) )
    lat = math.degrees(rlat)
    lon = math.degrees(rlon)
    zen = math.degrees(rzen)
    return ( lat, lon, zen )

  def locate(self,ix,iy):
    """ Convert ix,iy coords to lat,lon,zen angles

    PARAMETERS
      ix int : Current pixel x-coordinate
      iy int : Current pixel y-coordinate 

    RETURNS
      lat flt : Latitude [deg N]
      lon flt : Longitude [deg E]
      zen flt : Zenith angle [deg]

    DESCRIPTION 
      Uses GeoCal parameters to convert x,y to azi,ele angles 
      then uses locang to convert azi,ele to lat,lon,zen
      If GeoCal has not been set, returns (np.nan,np.nan,np.nan).
      May also return np.nan from locang if x,y beyond edge of disk
    """

    if self.cal:
      ele = ( iy - self.y0 ) / self.alpha
      azi = ( ix - self.x0 ) / self.beta
      return self.locang(ele,azi)
    else:
      return ( np.nan, np.nan, np.nan )

  def satang(self,lat,lon):
    """ Convert lat,lon angles to ele,azi,zen angles

    PARAMETERS
      lat flt : Latitude [deg N]
      lon flt : Longitude [deg E]

    RETURNS
      ele flt : Elevation angle [deg]
      azi flt : Azimuth angle [deg]
      zen flt : Zenith angle [deg]

    DESCRIPTION
      The inverse of locang
      Alpha is elevation and beta angle of rotation about inclined axis measured
      from the central vertical. Assumes spherical earth.
      Inputs/outputs in degrees, although internally converted to radians
    """

    # Convert lat,lon from degrees to radians
    rlat  = math.radians(lat)
    rlon  = math.radians(lon)
    # Height [km] of pixel above horizontal
    h2    = self.REARTH * math.sin(rlat)     
    # Distance [km] from earth's vertical axis
    r2    = self.REARTH * math.cos(rlat)     
    # Horizontal distance of pixel from satellite
    d3    = math.sqrt ( self.DIST**2 + r2**2 - 
                        2 * self.DIST * r2 * math.cos(rlon) )
    delta = math.atan ( h2 / d3 ) 
    gamma = math.asin ( r2 * math.sin(rlon) / d3 )
    rele  = math.atan ( math.tan(delta) / math.cos(gamma) )
    razi  = math.asin ( math.cos(delta) * math.sin(gamma) )
    rzen  = math.acos ( math.cos(razi)  * math.cos(rele) ) + \
            math.acos ( math.cos(rlat)  * math.cos(rlon) )         
    ele = math.degrees(rele)
    azi = math.degrees(razi)
    zen = math.degrees(rzen)
    return ( ele, azi, zen )
