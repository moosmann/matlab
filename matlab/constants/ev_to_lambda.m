function lambda = ev_to_lambda(eV)
% Converts energy given in electron volt (keV) into wavelength in metres.

% h = PlanckConstant('Js'); % J * s
% h_eV_s = 4.135667516e-15; eV * s
% c = speedOfLight; % m / s
%q      = 1.60217733e-16; %  J / eV [= h_J_s / h_eV_s ]
%lambda = h*c./(keV*q);

lambda = PlanckConstant('eVs') * SpeedOfLight / eV ;