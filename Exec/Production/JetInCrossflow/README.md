Case Description
----------------

Perform a {DNS, **LES} of a {nonreacting, **reacting} jet in crossflow where the jet issues into a {ducted, *channel} flow of vitiated combustion products at [atmospheric pressure]. The jet has a mole-basis composition of [70% H2, 18% N2 and 12% HE] at [300 K], a bulk inlet velocity of [42.2 m/s] with a {uniform profile, **uniform mean profile with turbulent fluctuations at an intensity of 5%, ****fully-developed turbulent structure}, and a diameter (d_j) of [3.175 mm]. The domain has a cross section of [24d_j] in the spanwise direction by [40d_j] in the streamwise direction and we are interested in the region from [30d_j] upstream of the jet to [50d_j] downstream of the jet. The walls are assumed to be {adiabatic, ***isothermal}. The cross flow has a bulk velocity of {19.1 m/s} with a {uniform profile, **uniform mean profile with turbulent fluctuations at an intensity of 5%, ****fully-developed turbulent structure} and a mole-basis composition of [12.91% O2, 76.11% N2, 3.66% CO2, 7.32% H2O and 0.0019629% OH] at [1236 K].

[] - Numerical Parameters that can be directly changed in the input file.
     Changing any of these should be the easiest level of change.
{} - Categorical options that may be changed with varying degress of difficulty. In all cases, the existing
     setup is the first option, and * are used to indicate the difficulty of alternative options:
     *: Easiest: simply change the value of an existing input parameters
     **: Somewhat harder: must add several input parameters not included in the current input file
     ***: Moderate difficulty: Slight changes needed to case source files
     ****: Extreme difficulty: Requires setting up and running a precursor simulation