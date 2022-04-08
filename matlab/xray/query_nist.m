function [energy, mass_att_coeff, energy_abs_coeff] = query_nist( element )

if nargin < 1
    element = 1;
end

% Get atomtic number
c = nist_material_constants_elements;
if ischar( element )
    element = strip( element );
    if numel( element ) <= 2
        m = strcmp(c{2}, element );
    else
        m = strcmpi(c{3}, element );        
    end
atomic_number = c{1}(m);
elseif isnumeric( element )
    atomic_number = element;
else
    error( '''element'' is neither an appropriate name, symbol or atomic number.' )
end


% Make URL
url = sprintf( 'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z%02u.html', atomic_number );

% Query webpage from NIST
o = weboptions('CertificateFilename','');
data = webread( url, o );

%% Parse webpage
num_char = numel( data );
% Header Line to skip
n = 0;
% current position in the character vector
pos = 0;
% allocation
energy = [];
mass_att_coeff = [];
energy_abs_coeff = [];
edge = [];
edge_energy = [];

% 
% while pos < num_char
%     [c, pos] = textscan( data, '%f %f %f', 1, 'HeaderLines', n);
%     if ~isempty( c{1} )
%         %disp( c )
%         energy = cat( 1, energy, c{1} );
%         mass_att_coeff = cat( 1, mass_att_coeff, c{2} );
%         energy_abs_coeff = cat( 1, energy_abs_coeff, c{3} );
%         
% %         [ce] = textscan( data, '%s %f %f %f',1, 'HeaderLines', n);
% %         disp(ce{1})
% %         disp(ce{2})
% %         if isnumeric(ce{2})
% %             edge = cat(1, edge, ce{1});
% %             edge_energy = cat(1, edge_energy, ce{2});
% %         end    
%         
%     end
%     n = n + 1;
% end

while pos < num_char
    [c, pos] = textscan( data, '%f %f %f', 'HeaderLines', n);
    if ~isempty( c{1} )
        %disp( c )
        n = n + numel( c{1} ) + 1;
        energy = cat( 1, energy, c{1} );
        mass_att_coeff = cat( 1, mass_att_coeff, c{2} );
        energy_abs_coeff = cat( 1, energy_abs_coeff, c{3} );
        
        [ce] = textscan( data, '%s %f %f %f',1, 'HeaderLines', n);
        disp(ce{1})
        disp(ce{2})
        if isnumeric(ce{2})
            edge = cat(1, edge, ce{1});
            edge_energy = cat(1, edge_energy, ce{2});
        end
    else        
        n = n + 1;
    end

end

 % energy is in MeV
 % mass_att_coeff in cm^2/g
 % energy_abs_coeff in cm^2/g
 
 % SI units conversion and uniqueness
 energy = ( energy ) * 1e6 ; % eV
 mass_att_coeff = ( mass_att_coeff ) / 10; % m^2 / kg
 energy_abs_coeff = ( energy_abs_coeff ) / 10; % m^2 / kg
 
 if numel( energy ) ~= numel( mass_att_coeff ) || numel( energy ) ~= numel( energy_abs_coeff )
     error( 'Mismatch of extracted arrays' )
 end
 