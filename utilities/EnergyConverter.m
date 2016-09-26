function lambda=EnergyConverter(keV)
% Converts Energy given in keV into wavelength in metres.    


%% MAIN
h      = 6.62606896e-34; % J * s
% h_J_s      = 6.62606896e-34;  
% h_eV_s = 4.135667516e-15; eV * s
c      = 299792458; % m / s
q      = 1.60217733e-16; %  J / eV [= h_J_s / h_eV_s ]
lambda = h*c./(keV*q);









