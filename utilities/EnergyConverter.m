function lambda=EnergyConverter(keV)
% Converts energy given in kilo electron volt (keV) into wavelength in
% metres.

% h = PlanckConstant('Js'); % J * s
% h_eV_s = 4.135667516e-15; eV * s
% c = speedOfLight; % m / s
%q      = 1.60217733e-16; %  J / eV [= h_J_s / h_eV_s ]
%lambda = h*c./(keV*q);

lambda = PlanckConstant('eVs') * speedOfLight / keV / 1000;









