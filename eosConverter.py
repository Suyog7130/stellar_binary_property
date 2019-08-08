
########################################################
###            Converting the EoS table              ###
########################################################

"""
04/08/2019:
This program is to get the EoS table in the desired format,

         * Column 1: Pressure (in MeV/fm^3)
         * Column 2: Energy Density (in MeV/fm^3)
         * Column 3: Number Desity (NA prob. in fm^{-3}, not req)

The 'eosaapr.dat' file by default has this format, whereas,
the 'eossly_original.dat' file has to changed to this one.
'eossly_original.dat' has original format of, 
Col1: Number Density (in /fm^3), Col2: Mass Density (in g/cm^3) and
Col3: Pressure (in MeV/fm^3).

"""

import math
import numpy as np
from numpy import pi
import scipy.constants as const


M_sun = 1.98847*10**30
G = const.G
c = const.c
k = G/c**2                #--this value is in N/kg^2.s^2 ~ m/kg
fact1 = G/c**4
log = math.log   


def main():
    
    L = np.loadtxt('eossly_original.dat')
    
    energy = L[:,1]     #--in g/cm^3.
    pressure = L[:,2]   #--in MeV/fm^3.
    num = L[:,0]
    
    energy *= 1000   #--in SI, kg/m^3.
    energy *= c**2   #--energy density, in SI, kg/m.s^2=Pa=J/m^3.
    
    conv1 = 1/(1.60218*10**32)
    energy *= conv1   #--in MeV/fm^3 (see Notes, Sec. 1)
    print(energy)
    print(len(energy), type(energy))
    
    #pressure *= 10   #--in SI, N/m^2=Pa.
    #pressure *= conv1   #--in MeV/fm^3 (see Notes, Sec. 1)
    #print(pressure)

    with open('eossly.dat', 'w') as f:
        for A, B, C  in zip(pressure, energy, num):
            f.write('{}\t{}\t{}\n'.format(A,B,C))


main()
        

   
#####################--End of Program--##################################
#########################################################################

