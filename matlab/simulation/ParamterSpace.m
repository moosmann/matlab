energy_min = 10;% keV
energy_max = 30;% keV
dist_min   = 0.01;% m
dist_max   = 5.00;% m
delta_min  = 1e-7;
delta_max  = 1e-6;
pixelsize  = 0.7e-6;% m


% energy/eV         delta             beta
% PMMA-formvar: chem. formula C5H8O2, density 1.19 gm/cm^3
%   1000.  0.000276933046  3.29402465E-05
%  10000.  2.67150153E-06  3.71373554E-09
%  20000.  6.66553376E-07  2.85367757E-10
%  30000.  2.96114848E-07  1.02135383E-10
%
energy = [energy_min,energy_max];
dist   = [dist_min,dist_max];
delta  = [delta_min,delta_max];

% lambda * distance
ld = dist'*EnergyConverter(energy);
ld = sort(ld(:));
ld
