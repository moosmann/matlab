function E = lambda_to_E(wavelength_in_metre)
% Convert energy in keV to wavelength in metre. The wavelength λ, in pm, can 
%be derived from the tabulated energy E, in keV, by the relationship λ = 1239.81/E.

% h      = 6.62606896e-34;
% c      = 299792458;
% q      = 1.60217733e-16;
% lambda = h*c./(keV*q);

E = 1.23984122e-09/wavelength_in_metre;