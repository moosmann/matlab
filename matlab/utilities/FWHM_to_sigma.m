function sigma = FWHM_to_sigma( FWHM )
% Convert FWHM of a Gaussian function to its standard deviaton sigma.

sigma = FWHM / sigma_to_FWHM( 1 );