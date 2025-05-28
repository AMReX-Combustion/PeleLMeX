## Single Droplet Evaporation Validation
This case compares results from PeleLMeX against experimental literature data.

| Test Case Name | $T_g$ [K] | $p_g$ [bar] | $T_d$ [K] | $d_d$ [$\mu m$] | Source |
| -------------- | --------- | ----------- | --------- | --------------- | ------ |
| Nomura         | 471, 741  | 1           | 298       | 700             | [1](https://doi.org/10.1016/S0082-0784(96)80344-4) |
| WongLin        | 1000      | 1.01325     | 315       | 2000            | [2](https://doi.org/10.1017/S0022112092003574)

1. H. Nomura, Y. Ujiie, H. J. Rath, J. Sato, and M. Kono, “Experimental study on high-pressure droplet evaporation using microgravity conditions,” Symposium (International) on Combustion, vol. 26, no. 1, pp. 1267–1273, Jan. 1996, doi: [10.1016/S0082-0784(96)80344-4](https://doi.org/10.1016/S0082-0784(96)80344-4).
2. S.-C. Wong and A.-C. Lin, “Internal temperature distributions of droplets vaporizing in high-temperature convective flows,” J. Fluid Mech., vol. 237, pp. 671–687, Apr. 1992, doi: [10.1017/S0022112092003574](https://doi.org/10.1017/S0022112092003574).

### Running Validation Cases
Compile the PeleLMeX case:
~~~
make TPL && make -j
~~~

The input file `input_general.inp` is set up for the `WongLin` case.  All other cases can be run by building different `case` objects in the `Validate.py` script.  Options include:
* `case = Nomura(471)`, `case = Nomura(741)`
* `case = WongLin()` 
