energy = 5 * 1e9 ; % eV

% Lorentz factor
gamma = energy / (electron_rest_mass('eV_c2') );

fprintf( '\n%10s:', 'E / GeV')
fprintf( '%10.2g', energy / 1e9) 

fprintf( '\n%10s:', 'gamma')
fprintf( '%10.2g', gamma) 

fprintf('\n')