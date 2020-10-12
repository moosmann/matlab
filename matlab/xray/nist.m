function [transmission, mac_tab, eac_tab, energy_tab, density, mat_const] = nist( element, energy_eV, thickness_m, show_plot )
% Return element specific interpolated transmission, the tabulated mass
% attenuation coefficients and tabulated energy absorption coefficient with
% according tabulated energies, the density and further material constants
% as given by NIST.
%
% ARGUMENTS:
% element : string or integer, name, symbol or atomic number
% energy_eV : vector, range of energies in eV
% thickness_m : material thickness in meter
% 
% RETURNS:
% transmission: vector of transmission values for given energy range
% mac_tab: vector of tabulated NIST mass attanuation coefficients μ/ρ in
%   m^2 / kg 
% eac_tab: vector of tabulated energy absorption coefficient  μ_en/ρ in m^2 /
%   kg 
% energy_tab: vector of tabulated energy values in eV where above NIST
%   values are given 
% density: scalar, density in  kg / m^2
% mat_const: 1 x 6 cell array of cell array of material parameters: { {atomic
%   number Z}, {symbol}, {name}, {Z/(mass A)}, {mean excitation energy I /
%   (eV)}, {density / (g/cm3)}}
% show_plot: bool, show transmission plot
%
% Written by J. Moosmann

%%  Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    element = 'hydrogen';
end
if nargin < 2
    energy_eV = (1:1:100)*1e3; % eV
end
if nargin < 3
    thickness_m = [10 100 1000] * 1e-3; % m
end
if nargin < 4
    show_plot = 0;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = nist_material_constants_elements;
num_thick = numel( thickness_m );
if ischar( element )
    element = strip( element );
    if numel( element ) <= 2
        m = strcmp(c{2}, element );
    else
        m = strcmpi(c{3}, element );        
    end    
else
    if isnumeric( element )
        m = element == c{1};
    else
        error( 'element is neither string or integer')
    end
end
mat_const = {c{1}(m), c{2}(m), c{2}(m), c{3}(m), c{4}(m), c{5}(m), c{6}(m)};
atomic_number = c{1}(m);

% X-ray mass attenuation coefficient mu/rho in m^2 / kg
name = c{3}(m);
name = lower( name{1} );
[energy_tab, mac_tab, eac_tab] = read_nist_txt( name );
if isempty( energy_tab )
    [energy_tab, mac_tab, eac_tab] = query_nist( atomic_number );
end

% density in kg / m^3
density = 1000 * c{6}(m);

if min( energy_eV ) < min( energy_tab )
    warning( 'Queried energy, %g eV, is below NIST tabulated values i.e. in [%g %g] eV', min( energy_eV ), min( energy_tab ), max( energy_tab ) )
end

if max( energy_eV ) > max( energy_tab )
    warning( 'Queried energy, %g eV, is below NIST tabulated values i.e. in [%g %g] eV', max( energy_eV ), min( energy_tab ), max( energy_tab ) )
end

% transmission I/I_0 = exp( - mu * t )
l = cell( [1, num_thick] );
for nn = num_thick:-1:1
    thick = thickness_m(nn);
    transmission(nn,:) = exp( - interp1( energy_tab, mac_tab, energy_eV ) .* density .* thick  );
    l{nn} = sprintf('thickness: %g mm',thick * 1000 );    
end

if show_plot
    figure( 'Name', sprintf( '%s', name ) )
    plot( energy_eV / 1000, transmission )
    xlabel('E / keV')
    ylabel('transmission')
    %l = sprintf('%u: %s, thickness: %g mm, density: %g g/cm^3', atomic_number, name, thickness_m * 1000, density / 1000 );
    t = sprintf('%u: %s, density: %g g/cm^3', atomic_number, name, density / 1000 );
    title( t )    
    legend( l )
    grid on
end
