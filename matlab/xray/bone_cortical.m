function [out, mac] = bone_cortical( energy, thickness_m )
% Return material constants for cortical bone
%
% Writen by J. Moosmann

% Bone, Cortical (ICRU-44)    
% <Z/A>	    0.51478	    
% I/(eV)  112.0   
% Composition/(Z: fraction by weight) 1: 0.034000 6: 0.155000 7: 0.042000 8: 0.435000 11: 0.001000 12: 0.002000 15: 0.103000 16: 0.003000 20: 0.225000

if nargin < 1
    energy = [20 30 40] * 1e3;
end
if nargin < 2
    thickness_m = 2e-3;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.density.value = 1.74e3; % Mg
out.density.unit = 'kg / m^3';

% elemental mass attenuation coefficients
[energy_bone_cortical, mass_att_coeff_bone_cortical] = read_nist_txt( 'bone_cortical' );

energy_min = min( energy_bone_cortical );
energy_max = max( energy_bone_cortical );
if min( energy ) < energy_min
    error( 'lowest queried energy, %g eV, must be above %g eV!', min( energy ), energy_min )
end
if max( energy ) > energy_max
    error( 'highest queried energy, %g eV, must be below %g eV!', max( energy ), energy_max )
end

% Compound coefficient obtained by additivity with mass weight wi
out.energy.value = energy;
out.energy.unit = 'eV';
out.mass_att_coeff.value = interp1( energy_bone_cortical, mass_att_coeff_bone_cortical, energy );
out.mass_att_coeff.unit = 'm^2 / kg';
mac = out.mass_att_coeff.value;

out.transmission.value = exp( - out.mass_att_coeff.value * out.density.value * thickness_m );
'm^2 / kg * kg / m^3 * m';
out.transmission.unit = 1;

out.absorption.value = 1 - out.transmission.value;
out.absorption.unit = 1;

%fprintf( '\n' )










