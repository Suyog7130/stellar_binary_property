# Stellar Binary Properties

This repository contains three python routines, ```tidal_precession.py```, ```tovsolver_onestar.py``` & ```tovsolver_k2.py```,  and two Equations of States (EoS), APR & SLY.

```tovsolver_onestar.py``` solves the TOV equation for a single star, going from its center to the surface. The radial variation of the mass can be visualised and tabulated. Inputs are the center energy density of the star in MeV/fm^3 and the EoS file. 

```tovsolver.py``` solves TOV equation for a given number of stars and the gets the Mass-Radius Curve for the chosen EoS. Additionally, it gets the value of the apsidal constant k2 for the Neutron Stars.

```tidal_precession.py``` calculates the value of the tidal component to the periastron precession of a binary. Input required are the two masses, radius of the first body and the apsidal constant.

A sample list of input values is, 
- Center energy density: 153.27 MeV/fm^3 (APR), 200.0 MeV/fm^3 (SLY).
- Number of Stars: 150 (APR), 110 (SLY).

The Equations of State, 'eosapr.dat' and 'eossly.dat' are saved in a three column format with the quatities in the order:
``` 
Pressure (MeV/fm^3) | Energy Density (MeV/fm^3) | Number Density (/fm^3) 
```

Our code can be used for other EoSs too. However, some EoSs might be available in the literature in different units and columns might be in a different order. Care should be taken to first convert the EoS into the above-mentioned form. As an example, eossly_original.dat is the original Sly EoS obtained from [SLY](https://github.com/thomascarreau/TOVsolver/tree/master/eos).

This has different columns as: ``` Number Density (in /fm^3) | Mass Density (in g/cm^3) | Pressure (in MeV/fm^3)  ```

The code ```eosConverter.py``` is used to convert it to the form required by our code ('eossly.dat'). 

TOV equation refers to the Tolman-Oppenheimer-Volkov Equation, which is a general relativistic hydrostatic equilibrium equation for a static, sperically symmetric star. For more look [here](https://www.wikiwand.com/en/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation).

See [Suyog Garg and Manjari Bagchi 2019 *Res. Notes AAS* **3** 125](https://doi.org/10.3847/2515-5172/ab3faf) for more details on the project.

[![DOI](https://zenodo.org/badge/197245843.svg)](https://zenodo.org/badge/latestdoi/197245843)
