% formulas: https://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html
% values: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/bone.html
% constants: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html

%Material 	  	           	   
% Bone, Cortical (ICRU-44)          
% <Z/A> 0.51478	
% I/(eV) 112.0 
% Density/(g/cm3) 1.920E+00	     1: 0.034000 6: 0.155000 7: 0.042000 8: 0.435000 11: 0.001000 12: 0.002000 15: 0.103000 16: 0.003000 20: 0.225000
% Composition/(Z: fraction by weight)

clear all
eV_to_J = 1 / 6.24e18;

E = (30:5:45)*1e3;
  
% Energy in eV
% X-ray mass attenuation coefficient mu/rho in m^2 / kg
[energy_bc, ~, mac_bc] = read_nist_txt( 'bone_cortical' );
[energy_wl, ~, mac_wl] = read_nist_txt( 'water_liquid' );

% thickness_bone in m
thickness_bone = ((1:3:4)*1e-3)'; 

% water density in kg / m^3
water_density = 1.0 * 1000;
thickness_water = ( 6e-3 - thickness_bone) / 2;
transmission_water = exp( - interp1( energy_wl, mac_wl, E ) * water_density .* thickness_water );

fprintf( '\ntransmission_water : %f', transmission_water )

% bone density in kg / m^3
bone_density = 1.920 * 1000; 

% transmission I/I_0 = exp( - mu * t )
mac_bc_cur = interp1( energy_bc, mac_bc, E );
absorption = 1 - exp( - mac_bc_cur * bone_density .* thickness_bone );
fprintf( '\n thickness bone : %f mm', thickness_bone * 1000)
fprintf( '\n absorption bone: %f', absorption * 100 )

% Flux area
flux = 1095734443447;
flux = 1.86e12;
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

full_exposure_time = 1; % seconds
mass = bone_density * area .* thickness_bone;

fprintf( '\n mass : %f mg ', mass * 1e6)

dose = transmission_water .* absorption * flux * full_exposure_time .* E * eV_to_J ./ mass;

fprintf( '\n dose : %f Gy', dose  )
figure('Name', 'Dose vs energy' )
plot( E, dose )
xlabel( 'energy' )
ylabel( 'dose' )
t = thickness_water;
l = {};
for nn=1:numel( t )
    l{nn} = sprintf( '%g', t(nn) );
end
legend( l )
axis tight
grid on

% figure( 'Name',  'log log plot: mass attenuation coefficient vs energy')
% loglog( E, mac_bc_cur *10 )
% xlabel( 'energy' )
% ylabel( 'mass attenuation coefficient' )
[energy_pet, ~, mac_pet] = read_nist_txt( 'pet' );
% thickness_bone in m
thickness_pet = 1e-3;
%  density in kg / m^3
pet_density = 1.3 * 1000;
transmission_pet = exp( - interp1( energy_pet, mac_pet, 45000 ) * pet_density .* thickness_pet );
fprintf( '\n transmission_pet : %f', transmission_pet )

fprintf( '\n')

mean_dose = 57.8534;
dose_level = [0.1 1 2 5] * 1e3;
dose_level = [0.1 1 2 5] * 1e3;

scan_time = ( dose_level / mean_dose );

num_samples = 3;
num_steps = 3;
total_scan_time =  num_samples * num_steps * cumsum( scan_time );

fprintf( '\n scan_times / s: ' )
fprintf( '\n  %f', round( scan_time ) )

fprintf( '\n accumulated scan_times / h: ' )
fprintf( '\n  %f', total_scan_time / 60 / 60 )



%fprintf( '', )


fprintf( '\n' )