
clear all

p(1).energy = 24.3e3; %eV
p(1).sample(1).label = 'pop18';
p(1).sample(1).folder = 'hnee18_pappel_tensionWood_001';

p(1).sample(2).label = 'pop18';


p(2).energy = 31e3; %eV
p(2).sample(1).label = 'pop20';
p(2).sample(2).label = 'pop21';
p(2).sample(3).label = 'pop22';
p(2).sample(4).label = 'pop23';


dir_fc = sprintf('/asap3/petra3/gpfs/p05/2018/data/11004450/processed/%s/flat_corrected', p(1).sample(1).folder);

exposure_time_scan = 12.5*60;% ~12 at 23keV, ~13 at 31keV;


flux = 1095734443447;
filename = '/asap3/petra3/gpfs/p05/2019/data/11005842/raw/dosetest/flux_window.tif';
im = double( imread( filename ) );
m = im > 1000;
pixelsize = 1.2801916e-6;
area = sum( m(:)) * pixelsize^2;
flux_density = flux / area;
fprintf( '\n flux : %g photons / s ', flux )
fprintf( '\n flux density : %g photons / s / m^2', flux_density )
fprintf( '\n pixelsize : %f micron', pixelsize * 1e6)
fprintf( '\n area: %f mm^2', area *1e6 )

thickness = 180e-6; %m
density = [0.35 0.5]' * g_per_ccm_to_kg_per_m(1); % kg/m^3
mass = density * area .* thickness;


fprintf( '\n mass : %f mg ', mass * 1e6)
energy = [23.4e3 31.0e3];
transmission = [0.75 0.8];
absorption = 0.01;
dose = transmission .* absorption * flux * exposure_time_scan .* energy * eV_to_J ./ mass;



fprintf('\n dose / kGy:\n' )

p = dose/1000;
p = cat(1, energy/1000, p);
p = cat(2, [0; density/1000], p);
disp(p)

fprintf('\n')
