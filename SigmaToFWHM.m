function FWHM = SigmaToFWHM(Sigma)
% Convert standard deviaton sigma of a Gaussian function to its FWHM.

FWHM = 2 * sqrt(2*log(2)) * Sigma;