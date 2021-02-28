# Weighted Spectral Difference

An R Package to Compare Molar Elipticity changes seen with Circular Dichroism across a temperature gradient using the Weighted Spectral Difference as introduced by Dinh et al., 2014.

* Dinh, N.N., Winn, B.C., Arthur, K.K. and Gabrielson, J.P., 2014. Quantitative spectral comparison by weighted spectral difference for protein higher order structure confirmation. Analytical biochemistry, 464, pp.60-62.

This package takes wide format dataframe(s) as inputs in which column 1 is the wavelength measured and subsequent columns are molar elipticity values (functionality for different units in future versions) at each wavelength at a given temperature. Column names of data should include a numeric value of temperature at which readings were taken (e.g. "ME-25C").

Seperate Dataframes should be used for seperate pH conditions and for seperate proteins (as applicable)

All dataframes can be input along with basic details of the data and calculated WSD values can be easily retrieved as a tidy data.frame suitable for plotting and further analysis.
