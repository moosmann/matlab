function sigma = FWHM_to_sigma( FWHM )
% Convert FWHM of a Gaussian function to its standard deviaton sigma.

sigma = FWHM / 2 * sqrt(2*log(2));