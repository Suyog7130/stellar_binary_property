
########################################################
###            Evaluating MR_curve & yR              ###
###                    SI Unit                       ###
########################################################

"""
This program can be used to plot the variation in Mass or 
the apsidal constant k2 with respect to Radius.

Runge-Kutta 4 is used to simultaneously solve equation 1, 2 
and 4 of README. 

The program can be run using the linux terminal. Python libraries:
'numpy', 'scipy.constants', 'scipy', 'matplotlib.pyplot'
and 'matplotlib.ticker' are required to use the program.

The example input values are:

        e0 = 153.27 (energy density in MeV/fm^3)
        N = 150      (no. of stars)

After this wait for some 30 seconds and then follow the instructions.
        
"""

import math
import numpy as np
from numpy import pi
import scipy.constants as const
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


M_sun = 1.98847*10**30
G = const.G
c = const.c
k = G/c**2                #--this value is in N/kg^2.s^2 ~ m/kg
fact1 = G/c**4
log = math.log   


#--Function for float range--#
def frange (start, stop, step):
    i = start
    while i<stop:
        yield i
        i += step
        
#--Plotting--#
def plotting (X, Y, marker, xname, yname, title):    
    fig, ax = plt.subplots(1, 1)
    ax.plot(X, Y, marker)
    ax.set_xlabel(xname, fontsize=15)
    ax.set_ylabel(yname, fontsize=15)
    ax.set_title(title, fontsize=15)
    ax.grid(True, which='major')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', which='major', labelsize=15, length=10)
    ax.tick_params(axis='both', which='minor', length=5)
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
#--------------------------------------------------------------------#   

#--Defining the first ODE (dP/dr) function--#
def f (x, y, z, s, args):
    k, c, e, dPdE = args
    A = z*c**2 + 4*pi*y*x**3
    num = -fact1 * (e+y) * A
    A = num / ( x**2 - 2*k*z*x )      #--A will have units J/m^4.
    return(A)

#--Defining the second ODE (dM/dr) function--#
def g (x, y, z, s, args):
    k, c, e, dPdE = args
    B = 4*pi*e*x**2/c**2              #--B will have units kg/m.
    return(B)  
    

#--Defining the third ODE (dyR/dr) function--#
def H (x, y, z, s, args):
    k, c, e, dPdE = args
    
    """#--conversions (refer the notes PDF, sec-1)--#
    x /= 1000                       #--in Km
    c /= 1000                       #--in Km/s
    k *= 1.476981739                #--in Km/M_sun
    z /= 1.98892*10**30             #--in M_sun 
    e *= 5.0278396266*10**(-28)
    y *= 5.0278396266*10**(-28)     #--in M_sun/km.s^2"""

    Fr = (x - 4*pi*(e-y)*(fact1)*x**3) / (x - 2*k*z)
    
    num1 = 4*pi*x * ( ( 5*e + 9*y + ((e+y)/dPdE) )*(fact1) - (6/(4*pi*x**2)) )
    den1 = x - 2*k*z
    part1 = num1/den1
    
    num2 = z*c**2 + 4*pi*y*x**3
    den2 = x**2 - 2*k*z*x
    part2 = 4*(fact1*num2/den2)**2
    Qr = part1 - part2
    
    H = (-1/x) * (s**2 + s*Fr + Qr*x**2) 
    
    return(H)

#--k2--#
def k2_value (yR, C):
    a = (8.0*C**5.0/5.0) * (1.0-2.0*C)**2.0 * (2.0+2.0*C*(yR-1.0)-yR)
    ba = 2.0*C * (6.0 - 3.0*yR + 3.0*C*(5.0*yR-8.0)) + 4.0*C**3.0 * (13.0 - 11.0*yR + C*(3.0*yR-2.0) + 2.0*C**2.0 * (1.0+yR))
    bb = 3.0*(1.0-2.0*C)**2.0 * (2.0 - yR + 2.0*C*(yR-1.0)) * log(1.0-2.0*C)
    return(a*(ba+bb)**(-1.0))

#---------------------------------------------------------------------------#    

#--Finding the matching value--#
def interpolateE (yi1, P, E): 

    for idx in range(1,len(P)-1):
        if yi1==P[idx]:
            e = E[idx]
            
        elif yi1>P[idx] and yi1<P[idx+1]:
            x0 = yi1
            x1, x2 = P[idx], P[idx+1]
            y1, y2 = E[idx], E[idx+1]
             
            e = ((y2-y1)*(x0-x1)/(x2-x1) ) + y1
            
    return(e)


def interpolateP (e0, P, E): 

    for idx in range(1,len(E)-1):
        if e0==E[idx]:
            p = P[idx]
            
        elif e0>E[idx] and e0<E[idx+1]:
            x0 = e0
            x1, x2 = E[idx], E[idx+1]
            y1, y2 = P[idx], P[idx+1]
             
            p = ((y2-y1)*(x0-x1)/(x2-x1) ) + y1
            
    return(p)
          
#------------------------------------------------------------------------#

#--RK4 Method--#
def rk4 (xi, yi, zi, si, h, args):

    k1 = f(xi, yi, zi, si, args)
    l1 = g(xi, yi, zi, si, args)
    j1 = H(xi, yi, zi, si, args)
    
    k2 = f(xi+h/2, yi+k1/2, zi+l1/2, si+j1/2, args)
    l2 = g(xi+h/2, yi+k1/2, zi+l1/2, si+j1/2, args)
    j2 = H(xi+h/2, yi+k1/2, zi+l1/2, si+j1/2, args)
    
    k3 = f(xi+h/2, yi+k2/2, zi+l2/2, si+j2/2, args)
    l3 = g(xi+h/2, yi+k2/2, zi+l2/2, si+j2/2, args)
    j3 = H(xi+h/2, yi+k2/2, zi+l2/2, si+j2/2, args)
    
    k4 = f(xi+h, yi+k3, zi+l3, si+j3, args)
    l4 = g(xi+h, yi+k3, zi+l3, si+j3, args)
    j4 = H(xi+h, yi+k3, zi+l3, si+j3, args)

    yi1 = yi + h/6*(k1 + 2*k2 + 2*k3 + k4)
    zi1 = zi + h/6*(l1 + 2*l2 + 2*l3 + l4)
    si1 = si + h/6*(j1 + 2*j2 + 2*j3 + j4)

    return(yi1, zi1, si1)
    

#--The TOV Solver Program--#
def eq_solve (L, e0, h):
    
    P = L[:, 0]                          #--first column is pressure
    E = L[:, 1]                          #--second column is energy density
    
    fact = 1.60218*10**32                
    P, E = fact*P, fact*E                #--converted to SI units MeV/fm^3 to Pa=J/m^2. 

    dP = np.gradient(P)
    dE = np.gradient(E)
    dPdE = dP/dE                         #--this is the array dP/de.

    xMax = 30000  
    x0 = 0.1      #--in meter.

    
    m0 = (4/3)*pi*x0**3*e0/c**2           #--this is in Kg.
    p0 = interpolateP(e0, P, E)           #--these are also in desired units. 
    dPdE0 = interpolateP(e0, dPdE, E)     #--this gives the starting value of dPdE.
    s0 = 2

    xi, yi, zi, si, ei, dPdEi = x0, p0, m0, s0, e0, dPdE0     #--initial conditions.
    
    for xi in frange(x0, xMax, h):
        args = [k, c, ei, dPdEi]
        yi, zi, si = rk4 (xi, yi, zi, si, h, args)
        
        if yi<0:
            break
        
        ei = interpolateE(yi, P, E)           #--thus e is already in the J/m^3.
        dPdEi = interpolateP(ei, dPdE, E)
    
    return(xi, yi, zi, si)

#---------------------------------------------------------------------------#

#--Main Program--#
def main ():

    try:
        EoS_file = input('\nPlease give the name of the EoS file (eosaapr.dat):\t')
        L = np.loadtxt(EoS_file)   #--Reading the EoS.
    except Exception:
        print('\n\tFILE NAME IS NOT CORRECT!!\n\n')
        raise

    #h = float(input('Give the step size:\t'))
    h = 10.0                         #--radial step size in meter.
    fact = 1.60218*10**32            #--to SI.
    
    e0 = float(input('Please give the value (eg: 153.27) of center energy density [MeV/fm^3]: '))
    print('\n e0 in J/m^2 is: ', e0*fact, '\n')
    
    N = int(input('Now give the number of stars (eg: 150): '))

    dE = 0
    radius, mass, yR = [], [], []
    
    for i in range(0, N):
        try:
            e0 = e0 + dE
            X, Y, Z, S = eq_solve(L, e0*fact, h)
            radius.append(X)
            mass.append(Z)
            yR.append(S)
            dE = 10.0
        except UnboundLocalError:
            print('\nThe max number of stars:', i, ' that can be plotted, is reached.')
            N = i
            break

    k2 = []
    for j in range(0, N):
        C = k*(mass[j]/radius[j])
        k2.append(k2_value(yR[j],C))
        
    radius, mass, yR, k2 = np.array(radius), np.array(mass), np.array(yR), np.array(k2)
    radius /= 1000
    mass /= M_sun  
    
    print('\n For the given EoS, we have:\n     Maximum Mass -> ',
          max(mass), ' M_sun', '\n	Corresponding Radius  -> ',
          float(radius[np.where(mass==max(mass))]), ' Km \n')   

    #--------------------------------------------------#
    ans1 = input('Do you want to store the output? (Y/N):\t')
    if ans1=='Y':
        fname1 = input('Give the table name in which to store the output:\t')
        with open(fname1, 'w') as f:
            for a, b, c  in zip(radius, mass, yR):
                f.write('{}\t{}\t{}\n'.format(a,b,c))
        
    ans2 = input('Do you want a plot? (Y/N):\t') 
    while ans2=='Y':
        query1 = input('On x-axis (radius/mass/k2):\t')
        query2 = input('On y-axis (radius/mass/k2):\t')
        query3 = input('Which marker (line/dot/dash):\t')
        title = query2 + '-' + query1 + ' Plot'
        
        if query1=='radius':
            x = radius
            xname = 'Radius $(Km)$'
        if query2=='radius':
            y = radius
            yname = 'Radius $(Km)$'
        if query1=='mass':
            x = mass
            xname = 'Mass $(M_{\odot})$'
        if query2=='mass':
            y = mass
            yname = 'Mass $(M_{\odot})$'
        if query1=='k2':
            x = k2
            xname = 'k2'
        if query2=='k2':
            y = k2
            yname = 'k2'
        if query3=='line':
            marker = 'c'
        if query3=='dot':
            marker = '.'
        if query3=='dash':
            marker = '--'    
            
        plotting(x, y, marker, xname, yname, title)
        ans2 = input('Do you want another plot? (Y/N):\t')    
    
    answer = str(input('Do you want to use another EoS (Y/N):\t'))
    while answer=='Y':
        answer='N'
        main()

    
main()
        

   
#####################--End of Program--##################################
#########################################################################
