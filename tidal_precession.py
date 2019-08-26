
#####################################################################
###      Evaluating the tidal component to the precession rate    ###
#####################################################################


import numpy as np
import pandas as pd
from numpy import pi
from scipy import integrate
import scipy.constants as const
from scipy import interpolate
import math
import matplotlib.pyplot as plt
from mpmath import *
from matplotlib.ticker import AutoMinorLocator
import pdfkit as pdf


M_sun = 1.98847*10**30 
R_sun = 6.957*10**8   
G = 6.67430*10**(-11)
c = 299792458             #--values taken from physics.nist.gov
k = G/c**2.0
fact = 1.476981739        #--to get G/c^2 in M_sun and Kms.
exp = math.exp
sin = math.sin
log = math.log            #--takes arguments ( , ) where the second one is the base, by default e.
    
#--------------------------------------------------------------------#   

def precessionEQ (Pb, Rns, aNS, Mbh, Mns, e, k2):
    #--the value found will be in SI, angle/s.
    fact1 = 30*pi/Pb
    fact2 = (Rns/aNS)**5
    fact3 = Mbh/Mns
    fact4 = (1 + 3/2*(e**2) + 1/8*(e**4)) / ((1-e**2)**5) 
    fact5 = k2
    return(fact1 * fact2 * fact3 * fact4 * fact5)

def aNSeq (Mbh, Mns, Pb):
    aR = ( (G*(Mbh+Mns)*(Pb**2)) / (4*(pi**2)) )**(1/3)
    aNS = ( Mbh/(Mbh+Mns) )*aR
    return(aR, aNS)



#--Main Program--#
def main ():
    
    print('\n Give the parameter values - ')
    Pb = float(input('\tPb (Hr):\t'))
    Rns = float(input('\tRns (km):\t'))
    Mbh = float(input('\tMbh (M_sun):\t'))
    Mns = float(input('\tMns (M_sun):\t'))
    e = float(input('\te:\t'))
    k2 = float(input('\tk2:\t'))
    
    #--Converting to SI--#
    Mbh, Mns = Mbh*M_sun, Mns*M_sun   #--masses in Kg.
    Rns = Rns*1000.0                  #--distances in meter.
    Pb = Pb*3600.0                    #--times in second.
        
    aR, aNS = aNSeq(Mbh, Mns, Pb)

    precession = precessionEQ(Pb, Rns, aNS, Mbh, Mns, e, k2)    
    precession = precession*1.807e+9   #--rad/sec to deg/year.
           
    print('We get - ')
    """print('Pb:\t', Pb, ' hours')
    print('Rns:\t', Rns/1000, ' km')
    print('Mbh:\t', Mbh/M_sun, ' M_sun')
    print('Mns:\t', Mns/M_sun, ' M_sun')
    print('e:\t', e)
    print('k2:\t', k2)
    print('aNS:\t', aNS/1000, ' km')"""
    print(' omegadot_tidal:\t', precession, ' deg/yr')
    
main()
        

    

  

############################--End of Program--##########################
########################################################################


