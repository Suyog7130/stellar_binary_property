# tovsolver

This repository contains two python routines, 'tovsolver_onestar.py' & 'tovsolver.py',  and two Equations of States (EoS), APR & SLY.

'tovsolver_onestar.py' solves the TOV equation for a single star, going from its center to the surface. The radial variation of the mass can be visualised and tabulated. Inputs are the center energy density of the star in MeV/fm^3 and the EoS file. 

'tovsolver.py' solves TOV equation for a given number of stars and the gets the Mass-Radius Curve for the chosen EoS. Additionally, it gets the value of the apsidal constant k2 for the Neutron Stars. The maximum value of k2 and the corresponding mass and radius values can be used to calculate the periastron precession rate of the eccentric binary if other parameters are known.

A sample list of input values is, 
- Center energy density: 153.27 MeV/fm^3 (APR), 200.0 MeV/fm^3 (SLY).
- Number of Stars: 150 (APR), 110 (SLY).

The Equations of State, 'eosapr.dat' and 'eossly.dat' are saved in a three column format with the quatities in the order:
 <b                           Pressure (MeV/fm^3) | Energy Density (MeV/fm^3) | Number Density (/fm^3)   >
 
 

TOV equation refers to the Tolman-Oppenheimer-Volkov Equation, which is a general relativistic hydrostatic equilibrium equation for a static, sperically symmetric star. For more look [here](https://www.wikiwand.com/en/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation).

<!--- See Garg & Manjari 2019 for more details. --->
