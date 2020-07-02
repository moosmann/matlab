function [energy, mass_att_coeff, energy_abs_coeff] = read_nist_txt( material )
% Returns the energy, the mass attenuation coefficient μ/ρ, and the energy
% absorption coefficient μen/ρ where ρ is the material density.
%
% ARGUMENTS
% material : string. Element or compound/mixture. One of the 'nist_.txt'
% files.
% 
% RETURNS
% energy : 1D vector, in eV
% mass_att_coeff : mass attenuation coefficient μ/ρ in m^2 / kg
% energy_abs_coeff : energye absorption coefficient  μen/ρ in m^2 / kg
%
% Written J. Moosman. First version. 2019-05-08

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    material = 'hydrogen';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 parpath = fileparts( mfilename( 'fullpath' ) );
 filename = sprintf( '%s/nist_%s.txt', parpath, material );
 fid = fopen( filename );
 c = textscan( fid, '%f %f %f', 'Delimiter', {'\n', '\r'}, 'CommentStyle', '%' );
 fclose( fid );
 
 % Energy         μ/ρ       μen/ρ 
% (MeV)       (cm2/g)    (cm2/g)
 
 % energy is in MeV
 % mass_att_coeff in cm^2/g
 % energy_abs_coeff in cm^2/g
 
 energy = c{1} * 1e6; % eV
 mass_att_coeff = c{2} / 10; % m^2 / kg
 energy_abs_coeff = c{3} / 10; % m^2 / kg
 