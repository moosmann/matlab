function fwhm = RMS_to_FWHM(rms)
% Assuming zero mean the root mean square of signal is identical to the
% standard deviation of the signal.
    fwhm = sigma_to_FWHM(rms);
