## Single Droplet Evaporation Validation
This case compares results from PeleLMeX against experimental literature data. Additional details and case descriptions are provided in the PelePhysics documentation at [https://amrex-combustion.github.io/PelePhysics/Spray.html#single-droplet-tests](https://amrex-combustion.github.io/PelePhysics/Spray.html#single-droplet-tests)

There are two general input files:
* `single-drop-evap-mp.inp`
* `single-drop-evap-gcm.inp`
Both are set up for the `WongLin` case with either the original PeleMP liquid properties or the GCM liquid properties.  The case with the PeleMP liquid properties can be run by compiling with `SPRAY_FUEL_NUM = 2` and
`SPRAY_GCM = FALSE`, then running:
~~~
mpirun -np 4 ./<PeleLMeX_EXE> single-drop-evap-mp.inp
~~~

Similarly, the case with the GCM liquid properties can be run by compiling with `SPRAY_FUEL_NUM = 2` and `SPRAY_GCM = TRUE`, then running:
~~~
mpirun -np 4 ./<PeleLMeX_EXE> single-drop-evap-gcm.inp
~~~

All cases provided in the PelePhysics documentation can be run by opening `Validate.py` and setting the `case_name` from the table in the PelePhysics documentation listed above, the model for the liquid properties, and the model for estimating saturated vapor pressure for the PeleMP model, which defaults to the Antoine fit. For example:
~~~
# Case to run
case_name = "Daif"

# Liquid properties model: "mp" or "gcm"
LiqPropsType = "mp"

# Psat model for PeleMP: "Antoine" or "Clasius-Clapeyron"
PeleMP_PsatModel = "Antoine"
~~~
then run
~~~
python Validate.py
~~~
Because many of the necessary options for each case must be set at compile time, the `Validate.py` script
will automatically recompile the code with the necessary options if a valid executable does not exist already.
The relevant compile time options are indicated in the executable name.

Case options include:
* `Nomura`
* `WongLin`
* `Daif`
* `RungeHep`, `RungeDec`, `RungeMix`, and `RungeJP8`

Users can compare the results from the various tests/configurations by running the `CompareLiqPropsType.py` script with the desired `case_name`. The script will search the current directory for all available data for each case.

Note that multicomponent evaporation is a work in progress as illustrated by the `RungeJP8` test case.

### Droplet Evaporation with Manifold-based Chemistry Models

This case setup also supports single droplet evaporation validation for spray modeling capability coupled
to manifold based chemistry models. These models replace the EOS and Transport property models from
PelePhysics with tabulated (or neural network) reduced-order representations. These models require additional
files containing the tabulated data and associated metadata, which are generated with the separate [CMLM
repository](https://github.com/NREL/cmlm). For example purposes, we include the necessary files to run the
`WongLin` case (`spray_wonglin.ctb` and `manifold_metadata_wonglin.text`). To run the sample case, which uses
PeleMP liquid properties with Antoine coefficients, first compile with `USE_MANIFOLD=TRUE`, `Manifold_Dim=1`,
`SPRAY_FUEL_NUM=1` and `SPRAY_GCM=FALSE`. Then run:
~~~
mpirun -np 4 ./<PeleLMeX_EXE> single-drop-evap-mp-manifold.inp
~~~

The `Validate.py` script can also be used to run any of the cases with manifold-based chemistry. This requires
an installed version of CMLM (including dependencies), which can be obtained within this directory using:
~~~
git clone git@github.com:NREL/cmlm.git
pip install -e cmlm
~~~
Then within `Validate.py` set `use_manifold = True` and the appropriate `cmlm_path` if downloaded elsewhere.
Note that the key script within CMLM that generates the tables for spray vaporization cases is located at
`cmlm/run_scripts/ctable/create_spray_table_nd.py`.