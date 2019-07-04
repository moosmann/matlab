function out = bone_cortical( energy, thickness_m )
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
out.density.value = 1.920e3;
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

out.transmission.value = exp( - out.mass_att_coeff.value * out.density.value * thickness_m );
'm^2 / kg * kg / m^3 * m';
out.transmission.unit = 1;

out.absorption.value = 1 - out.transmission.value;
out.absorption.unit = 1;


fprintf( '\n' )














% formulas: https://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html
% values: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/bone.html
% constants: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html

% Energy / MeV
E = [1.00000E-02 1.50000E-02 2.00000E-02 3.00000E-02 4.00000E-02 5.00000E-02 6.00000E-02 8.00000E-02];

% X-ray mass attenuation coefficient mu/rho / (cm^2/g)
mu_over_rho = [ 2.851E+01 9.032E+00 4.001E+00 1.331E+00 6.655E-01 4.242E-01 3.148E-01 2.229E-01];







