function [out, mac] = magnesium_gadolinium( energy_eV, mass_fraction_percent, thickness_m )
% Return material constants for a Mg-Gd alloy.
%
% Writen by J. Moosmann

% Arguments
if nargin < 1
    energy_eV = [10 30 40] * 1e3;
end
if nargin < 2
    mass_fraction_percent = 10; % mass percent
end
if nargin < 3
    thickness_m = 2e-3;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.formula = sprintf( 'Mg%uGd', mass_fraction_percent );
if mass_fraction_percent == 5
    out.density.value = 1811;
elseif mass_fraction_percent == 10
    out.density.value = 1876;
else
    error( 'Density of alloy ratio unknown. Extrapolate.' )
end

out.density.unit = 'kg / m^3';
out.density_Mg.value = 1738; % Wikipedia
out.density_Mg.unit = 'kg / m^3';
out.density_Gd.value = 7900; % Wikipedia
out.density_Gd.unit = 'kg / m^3';

% elemental mass attenuation coefficients
[energy_Mg, mass_att_coeff_Mg] = read_nist_txt( 'magnesium' );
[energy_Gd, mass_att_coeff_Gd] = read_nist_txt( 'gadolinium' );

% Check energy range
energy_min = min( [ energy_Mg' energy_Gd' ] );
energy_max = max( [ energy_Mg' energy_Gd' ] );
if min( energy_eV ) < energy_min
    error( 'lowest queried energy_eV, %g eV, must be above %g eV!', min( energy_eV ), energy_min )
end
if max( energy_eV ) > energy_max
    error( 'highest queried energy_eV, %g eV, must be below %g eV!', max( energy_eV ), energy_max )
end


% 10 refers to 10 mass-percent. With a standard atomic weight of 

% weight percentage
out.molar_mass_Mg.value = 24.305; % Wikipedia
out.molar_mass_Mg.unit = 'g / mol';
out.molar_mass_Gd.value = 157.3;
out.molar_mass_Gd.unit = 'g / mol';

% mass fraction: w_i = x_i M_i / M, x_i: mole fraction, M: molar mass
% average
w_Mg = ( 1 - mass_fraction_percent / 100 );
w_Gd = mass_fraction_percent / 100;

% Average molar mass of mixtures M
% M = sum_i x_i M_i, x_i: mole fraction, M_i: molar mass of componten i
% or M = 1 / ( sum_i w_i / M_i ), w_i: mass fraction
out.molar_mass.value = 1 / (w_Mg / out.molar_mass_Mg.value + w_Gd / out.molar_mass_Gd.value );
out.molar_mass.unit = 'g / mol';

% mole fraction
out.molar_fraction_Mg = w_Mg * out.molar_mass.value / out.molar_mass_Mg.value;
out.molar_fraction_Gd = w_Gd * out.molar_mass.value / out.molar_mass_Gd.value;

% Stochiometric ratio
% N_Gd * u_Gd = X / 100 * ( N_Mg * u_Mg + N_Gd * u_Gd )
out.stochiometric_ratio.value = out.molar_fraction_Mg / out.molar_fraction_Gd;
%out.stochiometric_ratio.value = out.molar_mass_Gd.value * (1 - X / 100) / (out.molar_mass_Mg.value * X / 100 );
out.stochiometric_ratio.unit = '1 = N_Mg / N_Gd';

% Compound coefficient obtained by additivity with mass weight wi
out.energy.value = energy_eV;
out.energy.unit = 'eV';
out.mass_att_coeff.value = ...
      w_Mg * interp1( energy_Mg, mass_att_coeff_Mg, energy_eV ) ...
    + w_Gd * interp1( energy_Gd, mass_att_coeff_Gd, energy_eV );
out.mass_att_coeff.unit = 'm^2 / kg';
mac = out.mass_att_coeff.value;

out.transmission.value = exp( - out.mass_att_coeff.value * out.density.value * thickness_m );
'm^2 / kg * kg / m^3 * m';
out.transmission.unit = 1;

out.absorption.value = 1 - out.transmission.value;
out.absorption.unit = 1;

%fprintf( '\nweights: Cd %f, W: %f, O: %f (total: %f)', w_Mg, w_Gd, w_O, w_Mg + w_Gd + w_O )
%fprintf( '\n' )
