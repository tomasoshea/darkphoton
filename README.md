# darkphoton

Respository to host code, data and plots for the calculation of new limits on the dark photon parameters that would be obtained by the future IAXO helioscope.

In its current form, the solar flux of transversely and longitudinally polarised dark photons are read from data/output_T.dat and data/output_L.dat respectively, calculated in Javier Redondo's [Atlas of Solar Hidden Photon Emission](https://arxiv.org/abs/1501.07292). The solar dark photon flux from plasmon mixing and nuclear sources are calculated directly in the code as specified below. The ability to directly calculate output_T.dat and output_L.dat exists in **darkphoton.h** but is not currently used.

Note that all the code in the DPlimits folder is taken from [cajohare/AxionLimits](https://github.com/cajohare/AxionLimits) and slightly modified to include the calculated IAXO limits. This is used to create the nice looking limits plots.

* **filter.py** calculates solar parameters from the solar model [AGSS09](https://arxiv.org/abs/0909.2668) and outputs these as .dat files for use in C++ integration

*  **darkphoton.h** contains useful code to be run in the .cpp files

*  **fluxAtlas.cpp** outputs the back-converted photon flux in the vacuum and gas runs for thermally produced transverse dark photons

*  **LMixing.cpp** does the same for transverse dark photons originating from longitudinal solar plasmons

*  **deexcitation.cpp** does the same for dark photons created in gamma transitions

*  **pp-chain.cpp** does the same for dark photons created in nuclear fusion and e+ e- annihilation

*  **fluxAltas-pureL.cpp** gives the theoretical flux were an IAXO-like helioscope sensitive to longitudinal oscillations in the buffer gas

* **Bfields.cpp** calculates the solar B-fields from a model taken from [New Solar Seismic Models and the Neutrino Puzzle](https://arxiv.org/abs/astro-ph/0203107v1)

* The remaining c++ codes output parameters that it is useful to plot, for example **Gamma.cpp** outputs the absorption length in the IAXO buffer gas as a function of pressure or photon energy

* The remaining python scripts are to plot the various quantities calculated in the C++ code 
