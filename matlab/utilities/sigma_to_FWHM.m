function FWHM = sigma_to_FWHM(sigma)
% Convert standard deviaton sigma of a Gaussian function to its FWHM.

FWHM = 2 * sqrt(2*log(2)) * sigma;