function out = PEEK( energy_eV, thickness_m )
% Return material constants for PEEK
%
% Written by J. Moosmann

if nargin < 1
    energy_eV = [20 30 40] * 1e3;
end
if nargin < 2
    thickness_m = 5e-3;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out.energy.value = energy_eV;
out.energy.unit = 'eV';
out.thickness.value = thickness_m;
out.thickness.unit = 'm';

out.formula = 'C19H12O3';
out.density.value = 1320;
out.density.unit = 'kg / m^3';
out.molar_mass.value = 288.30;
out.molar_mass.unit = 'g / mol';

% elemental mass attenuation coefficients
[energy_C, mass_att_coeff_C] = read_nist_txt( 'carbon' );
[energy_H, mass_att_coeff_H] = read_nist_txt( 'hydrogen' );
[energy_O, mass_att_coeff_O] = read_nist_txt( 'oxygen' );

energy_min = min( [energy_C' energy_H' energy_O'] );
energy_max = max( [energy_C' energy_H' energy_O'] );
if min( energy_eV ) < energy_min
    error( 'lowest queried energy_eV, %g eV, must be above %g eV!', min( energy_eV ), energy_min )
end
if max( energy_eV ) > energy_max
    error( 'highest queried energy_eV, %g eV, must be below %g eV!', max( energy_eV ), energy_max )
end

% weight percentage
molar_mass_C.value = 12.011;
molar_mass_C.unit = 'g / mol';
molar_mass_H.value = 1.008;
molar_mass_H.unit = 'g / mol';
molar_mass_O.value = 15.999;
molar_mass_O.unit = 'g / mol';
molar_mass_PEEK_calculated.value = 19 * molar_mass_C.value + 12 * molar_mass_H.value + 3 * molar_mass_O.value;
molar_mass_PEEK_calculated.unit = 'g / mol';
molar_mass_total = out.molar_mass.value;
fprintf( '\nmolar mass:\n  tabulated: %g %s\n calculated: %g %s', out.molar_mass.value, out.molar_mass.unit, molar_mass_PEEK_calculated.value, molar_mass_PEEK_calculated.unit ) 
w_C = 1 * molar_mass_C.value / molar_mass_total;
w_W = 1 * molar_mass_H.value / molar_mass_total;
w_O = 4 * molar_mass_O.value / molar_mass_total;

% Compound coefficient obtained by additivity with mass weight wi
out.mass_att_coeff.value = ...
      w_C * interp1( energy_C, mass_att_coeff_C, energy_eV ) ...
    + w_W * interp1( energy_H, mass_att_coeff_H, energy_eV ) ...
    + w_O * interp1( energy_O, mass_att_coeff_O, energy_eV );
out.mass_att_coeff.unit = 'm^2 / kg';

out.transmission.value = exp( - out.mass_att_coeff.value * out.density.value * thickness_m );
'm^2 / kg * kg / m^3 * m';
out.transmission.unit = 1;

out.absorption.value = 1 - out.transmission.value;
out.absorption.unit = 1;

%fprintf( '\nweights: Cd %f, W: %f, O: %f (total: %f)', w_Cd, w_W, w_O, w_Cd + w_W + w_O )
%fprintf( '\n' )