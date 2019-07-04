function out = cadmium_tungstate( energy, thickness_m )
% Return material constants for cadmium tungstate (CdWO4).
%
% Writen by J. Moosmann

if nargin < 1
    energy = [20 30 40] * 1e3;
end
if nargin < 2
    thickness_m = 50e-6;
end

out.formula = 'CdWO4';
out.density.value = 7.9e3;
out.density.unit = 'kg / m^3';
out.molar_mass.value = 360.25;
out.molar_mass.unit = 'g / mol';

% elemental mass attenuation coefficients
[energy_Cd, mass_att_coeff_Cd] = read_nist_txt( 'cadmium' );
[energy_W, mass_att_coeff_W] = read_nist_txt( 'tungsten' );
[energy_O, mass_att_coeff_O] = read_nist_txt( 'oxygen' );

energy_min = min( [energy_Cd' energy_W' energy_O'] );
energy_max = max( [energy_Cd' energy_W' energy_O'] );
if min( energy ) < energy_min
    error( 'lowest queried energy, %g eV, must be above %g eV!', min( energy ), energy_min )
end
if max( energy ) > energy_max
    error( 'highest queried energy, %g eV, must be below %g eV!', max( energy ), energy_max )
end

% weight percentage
molar_mass_Cd.value = 112.414;
molar_mass_Cd.unit = 'g / mol';
molar_mass_W.value = 183.84;
molar_mass_W.unit = 'g / mol';
molar_mass_O.value = 15.999;
molar_mass_O.unit = 'g / mol';
molar_mass_total = out.molar_mass.value;
w_Cd = 1 * molar_mass_Cd.value / molar_mass_total;
w_W = 1 * molar_mass_W.value / molar_mass_total;
w_O = 4 * molar_mass_O.value / molar_mass_total;



% Compound coefficient obtained by additivity with mass weight wi
out.energy.value = energy;
out.energy.unit = 'eV';
out.mass_att_coeff.value = ...
    w_Cd * interp1( energy_Cd, mass_att_coeff_Cd, energy ) ...
    + w_W * interp1( energy_W, mass_att_coeff_W, energy ) ...
    + w_O * interp1( energy_O, mass_att_coeff_O, energy );
out.mass_att_coeff.unit = 'm^2 / kg';

out.transmission.value = exp( - out.mass_att_coeff.value * out.density.value * thickness_m );
'm^2 / kg * kg / m^3 * m';
out.transmission.unit = 1;

out.absorption.value = 1 - out.transmission.value;
out.absorption.unit = 1;

%fprintf( '\nweights: Cd %f, W: %f, O: %f (total: %f)', w_Cd, w_W, w_O, w_Cd + w_W + w_O )
%fprintf( '\n' )
