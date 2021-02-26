'''
These functions will be used to convert between radiances and temperatures
'''
import math
import numpy as np

#------------------------------------------------------------------------
# LOCAL CONSTANTS
# These are fundamental constants used in the calculations
# Lifted straight from the Image class
H = 6.626e-34       # Planck constant       [m2 kg / s]
C = 2.9979e8         # Speed of light        [m / s]
K = 1.3806e-23       # Boltzmann constant    [m2 kg /s2 /K]
R1 = H * C / K     # Intermediate Constant [m K]             0.014413
R2 = 2 * H * C**2  # Intermediate Constant [m4 kg / s3]      1.19e-16

#------------------------------------------------------------------------
# RADIANCE FUNCTION
# This is Plank's function, used to calculate the radiance of a blackbody at a bgiven wavelength and temperature
'''
INPUTS:
    T [float] : temperature (kelvin) for the radiance to be caluclated
    wavelength [float] : wavelength at which Plank function should be evaluated (microns)

OUTPUTS:
    R [float] : calculated radiance (W/(m2 sr um))
'''
def B(T,wavelength):
    w = wavelength * 1.0e-6  # convert microns to metres
    R = R2 / w**5 / 1e6 / (np.exp(R1 / (T * w))-1)
    return R

#-------------------------------------------------------------------------
# TEMPERATURE FUNCTION
# Uses the inverse of Plank's function to calculate a temperature given a radiance
'''
INPUTS:
    R [float] : radiance for which we want the temperature (W/(m2 sr um))
    wavelength [float] : wavelength for distribution to be evaluated at (microns)

OUTPUTS:
    t [float] : calculated temperature (kelvin)
'''
def T(R,wavelength):
    w = wavelength * 1.0e-6  # convert microns to metres
    t = R1 / w / np.log( 1.0 + R2/(w**5 * R*1e6) )
    return t

#--------------------------------------------------------------------------
# BRIGHTNESS ERROR FUNCTION
# calculataes the error in a brightness value given the error in a temperature value
'''
INPUTS:
    T [float] : temperature value
    wavelength [float] : wavelength inmicrons
    dT [float] : error in temperature that will give error in radiance

OUTPUTS:
    err [float] : error in radiance
'''
def Berr(T,wavelength,dT):
    w = wavelength * 1.0e-6  # convert microns to metres
    b = B(T,wavelength)

    err = b*R1 / w / T**2 / (1-np.exp( -R1 / (w * T) )) * dT
    
    return err

#--------------------------------------------------------------------------
# TEMPERATURE ERROR FUNCTION
# calculataes the error in a temperature value given the error in a brightness value
'''
INPUTS:
    R [float] : temperature value
    wavelength [float] : wavelength in microns
    dR [float] : error in temperature that will give error in radiance

OUTPUTS:
    err [float] : error in radiance
'''
def Terr(R,wavelength,dR):
    w = wavelength * 1.0e-6  # convert microns to metres
    t = T(R,wavelength)

    err = t / R / ( 1 + (w**5 * R / R2) ) / np.log( 1.0 + R2/(w**5 * R) ) * dR
    
    return 3*err
    

