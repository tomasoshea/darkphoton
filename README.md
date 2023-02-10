# darkphoton

Respository to host code, data and plots for the calculation of new limits on the dark photon parameters that would be obtained by the future IAXO helioscope.

Note that all the code in the DPlimits folder is taken from [cajohare/AxionLimits](https://github.com/cajohare/AxionLimits) and slightly modified to include the calculated IAXO limits. This is used to create the nice looking limits plots.

* **filter.py** calculates solar parameters from the solar model [AGSS09](https://arxiv.org/abs/0909.2668) and outputs these as .dat files for use in C++ integration

* **IAXOstats.py** calculates the limit on the flux for a given data taking period, confidence level and background flux

*  **darkphoton.h** contains all the useful code to be run in the .cpp files

* **tPlasmon.cpp** and **tPlasmonGas.cpp** calculate the limits for transverse DPs in the IAXO vacuum and gas runs respectively

* **lMixingRes.cpp** and **lMixingResGas.cpp** do the same for the case of transverse DPs emitted from longitudinal plasmons

* **pureL.cpp** calculates the limits that could be obtained by the direct detection of longitudinal DPs

* **Espectrum.cpp** and its variants calculate the flux spectrum for each method of DP production

* **Bfields.cpp** calculates the solar B-fields from a model taken from [New Solar Seismic Models and the Neutrino Puzzle](https://arxiv.org/abs/astro-ph/0203107v1)

* The remaining python scripts are to plot the various quantities calculated in the C++ code
