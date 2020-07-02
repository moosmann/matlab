% Parameter settings for spiral CT
% Felix : welche camera, which binning, which speed at which binning, min
% exposure times
% Fabian : diode

eV_to_J = 1 / 6.24e18;

exp_time = 0;
bin_onchip = 2;
bin_raw = 3;
bin_total =  bin_onchip * bin_raw;
magn = 5;
% IBL
%flux = 5.4e10 / 65 / 48; % photons / s / mm2
% Tomcat 30 keV
energy = 30e3; % eV

% Dose
% Water
water_density = 1.0 * 1000; % water density in kg / m^3
[energy_wl, ~, mac_wl] = read_nist_txt( 'water_liquid' );
mac_water = interp1( energy_wl, mac_wl, energy );
% Bone
[energy_bc, ~, mac_bc] = read_nist_txt( 'bone_cortical' );
mac_bone = interp1( energy_bc, mac_bc, energy );
% Fluxe
flux = 2.1e11 * 1e6; % photons / s / mm2
entrance_doserate_water = flux * energy * eV_to_J * mac_water;
entrance_doserate_bone = flux * energy * eV_to_J * mac_bone;

% Camera info
p = pp_info;
cam = p.camera(3);
fps = cam.framerate_12bit_onchipbin2x_Hz;
exp_time_min = 1 / fps * 1000;
switch magn
    case 5
        fov = cam.fov_5x_mm;
    case 10
        fov = cam.fov_10x_mm;
    case 20
        fov = cam.fov_20x_mm;
    case 40
        fov = cam.fov_40x_mm;
end
fov = round( fov, 1);
shape = cam.pixels;
shape_bin = floor( shape / bin_total );
pix_hor = shape_bin(1);
pix_vert = shape_bin(2);
eff_pix_bin = 1 / magn * cam.pixel_size_micron * bin_total;
% Numbper of projection
num_proj = round( pix_hor * pi / 2 );
% Scan time
scan_time = num_proj * exp_time_min / 1000;
% Rotation blurring
blur_ang = pi * pix_hor /2 / num_proj ;
% Translation  blurring: During 180 degree the sample must move less than
% one detector height
blur_vert = pix_vert / num_proj;

fprintf( '\n camera : %s', cam.name )
fprintf( '\n pixel binned : %u x %u', shape_bin )
fprintf( '\n fov : %g mm x %g mm', fov )
fprintf( '\n effective pixel size binned : %f micron', eff_pix_bin )
fprintf( '\n num_proj (pix_hor * pi / 2) : %u', num_proj )
fprintf( '\n max angular blur : %f * eff_pix_bin', blur_ang )
fprintf( '\n vertical blurring : %f * eff_pix_bin', blur_vert )
fprintf( '\n fps : %u Hz', fps )
fprintf( '\n full frame exposure time : %f ms', exp_time_min )
fprintf( '\n scan time / half turn : %f s', scan_time )
fprintf( '\n flux density : %.2g photons / s / mm2', flux * 1e-6 )

fprintf( '\n entrance dose water : %.3g Gy / s', entrance_doserate_water )
fprintf( '\n entrance dose bone : %.3g Gy / s', entrance_doserate_bone )

fprintf( '\n\n' )