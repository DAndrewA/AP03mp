'''
This function is to caluclate the characteristic wavelength of each filter.
This will allow us to convert our radiances into brightness temperatures.
'''
from pathlib import Path

### INPUTS ###
# channel [int] : the channel the characteristic wavelength is being calculated for
# direc [Path] : the relative directory the files are kept in. This should be ./pectralResponses

### OUTPUTS ###
# wavelength [float] : the characteristic wavelength of the filter in microns 

def specWav(channel,direc):
    # opens text file for given spectral response
    filename = "ch" + str(channel) + "Response.txt"
    f = open(direc / filename,'r')
    contents = f.read()
    f.close()

    # splits file into each row of data
    contents = contents.split("\n")
    contents = contents[6:]

    spectralResponse = 0
    normalisation = 0
    for row in contents:
        try:
            data = row.split()
            spectralResponse += float(data[0])*float(data[1])
            normalisation += float(data[1])
        except:
            #print("Problem with line reading: " + str(row))
            a = 1

    wavelength = spectralResponse/normalisation
    return wavelength
