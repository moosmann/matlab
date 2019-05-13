function out = PEEK
% Return material constants for PEEK

out.formula = 'C19H12O3';
out.density.value = 1320;
out.density.unit = 'kg / m^3';

% elemental mass attenuation coefficients
[~, C] = read_nist_txt( 'carbon' );
[~, H] = read_nist_txt( 'hydrogen' );
[~,O] = read_nist_txt( 'oxygen' );

% Compound coefficient obtained by additivity with mass weight wi
out.mass_att_coeff.value = wC * C + wH * H + wO * O;
out.mass_att_coeff.unit = 'cm^2 / g';

%