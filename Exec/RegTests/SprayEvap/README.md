## Single Droplet Evaporation Validation
This case compares results from PeleLMeX against experimental literature data.

| Test Case Name | $\bf{T_g}$ [K] | $\bf{p_g}$ [bar] | $\bf{T_d}$ [K] | $\bf{d_d}$ [μm] | $\bf{\Delta u}$ [m/s] | Source |
| -------------- | -------------- | ---------------- | -------------- | --------------- | ------------------- | ------ |
| Nomura         | 471, 741  | 1           | 298       | 700   | 0.0   | [1](https://doi.org/10.1016/S0082-0784(96)80344-4) |
| WongLin        | 1000      | 1.01325     | 315       | 1961  | 0.385 | [2](https://doi.org/10.1017/S0022112092003574)
| Daif           | 348       | 1.01325     | 291.4     | 1334  | 3.1   | [3](https://doi.org/10.1016/S0894-1777(98)10035-3)


### Running Validation Cases
Compile the PeleLMeX case:
~~~
make TPL && make -j
~~~

The input file `input_general.inp` is set up for the `WongLin` case.  This can be run using:
~~~
mpirun -np 4 ./<PeleLMeX_EXE> input_general.inp
~~~

All other cases can be run by building different `case` objects in the `Validate.py` script, and running: 
~~~
python Validate.py
~~~

Case options include:
* `case = Nomura(471)`, `case = Nomura(741)`
* `case = WongLin()`
* `case = Daif()`

### References
1. H. Nomura, Y. Ujiie, H. J. Rath, J. Sato, and M. Kono, “Experimental study on high-pressure droplet evaporation using microgravity conditions,” Symposium (International) on Combustion, vol. 26, no. 1, pp. 1267–1273, Jan. 1996, doi: [10.1016/S0082-0784(96)80344-4](https://doi.org/10.1016/S0082-0784(96)80344-4).
2. S.-C. Wong and A.-C. Lin, “Internal temperature distributions of droplets vaporizing in high-temperature convective flows,” J. Fluid Mech., vol. 237, pp. 671–687, Apr. 1992, doi: [10.1017/S0022112092003574](https://doi.org/10.1017/S0022112092003574).
3. A. Daı̈f, M. Bouaziz, X. Chesneau, and A. Ali Chérif, “Comparison of multicomponent fuel droplet vaporization experiments in forced convection with the Sirignano model,” Experimental Thermal and Fluid Science, vol. 18, no. 4, pp. 282–290, Dec. 1998, doi: [10.1016/S0894-1777(98)10035-3](https://doi.org/10.1016/S0894-1777(98)10035-3).
