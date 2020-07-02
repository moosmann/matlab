% P05 paramter
clear all

%% Source PETRA III

% distance source sample
p.coherence.dist_source_dcm = 50.9; % m
p.coherence.dist_source_sample = 82.7; % m

% coherence length
p.coherence.sigma_h = 38e-6; % m
p.coherence.sigma_v = 6e-6; % m
p.coherence.FWHM_h = sigma_to_FWHM( p.coherence.sigma_h );
p.coherence.FWHM_v = sigma_to_FWHM( p.coherence.sigma_v );

%% Optics
p.opt.magn = [5 10 20 40];
p.opt.producer = 'Jena Optik';

%% CCD camera
p.camera.ccd.name = 'EHD SciCam SC09000';
p.camera.ccd.sensor = 'KAF-09000';
p.camera.ccd.dynamic = '16 bit';
p.camera.ccd.pixels = [3056 3056];
p.camera.ccd.pixel_size = 12e-6;
p.camera.ccd.effective_pixel_size_micron = p.camera.ccd.pixel_size * 1e6 ./ p.opt.magn;
p.camera.ccd.pixel_size_micron = p.camera.ccd.pixel_size * 1e6;
p.camera.ccd.fov = p.camera.ccd.pixel_size * p.camera.ccd.pixels(1) ./ p.opt.magn;
p.camera.ccd.fov_mm = p.camera.ccd.fov * 1e3;

%% CMOS KIT camera
p.camera.cmos(1).name = 'KIT';
p.camera.cmos(1).sensor = 'CMOSIS CMV20000';
p.camera.cmos(1).dynamic = '12 bit';
p.camera.cmos(1).framerate_at_12bit = 30;
p.camera.cmos(1).pixels = [5120 3840];
p.camera.cmos(1).pixel_size = 6.4e-6;
p.camera.cmos(1).effective_pixel_size_micron = p.camera.cmos(1).pixel_size * 1e6 ./ p.opt.magn;
p.camera.cmos(1).chip_size_mm = p.camera.cmos(1).pixels .* p.camera.cmos(1).pixel_size * 1e3;
p.camera.cmos(1).pixel_size_micron = p.camera.cmos(1).pixel_size * 1e6;
p.camera.cmos(1).fov = p.camera.cmos(1).pixel_size * p.camera.cmos(1).pixels(1) ./ p.opt.magn;
p.camera.cmos(1).fov_mm = p.camera.cmos(1).fov * 1e3;

%% CMOS KIT camera
p.camera.cmos(2).name = 'XIMEA';
p.camera.cmos(2).sensor = 'CMOSIS ';
p.camera.cmos(2).dynamic = '12 bit';
p.camera.cmos(2).framerate_at_12bit = 30;
p.camera.cmos(2).pixels = [5120 3840];
p.camera.cmos(2).pixel_size = 6.4e-6;
p.camera.cmos(2).effective_pixel_size_micron = p.camera.cmos(2).pixel_size * 1e6 ./ p.opt.magn;
p.camera.cmos(2).chip_size_mm = p.camera.cmos(2).pixels .* p.camera.cmos(2).pixel_size * 1e3;
p.camera.cmos(2).pixel_size_micron = p.camera.cmos(2).pixel_size * 1e6;
p.camera.cmos(2).fov = p.camera.cmos(2).pixel_size * p.camera.cmos(2).pixels(1) ./ p.opt.magn;
p.camera.cmos(2).fov_mm = p.camera.cmos(2).fov * 1e3;


%% Stages
p.stage.rotation.vendor = 'Aerotech';
p.stage.rotation.travel_range_z = [-10,10];

p.stage.base.vendor = 'Aerotech';
p.stage.base.type = 'tripod';
p.stage.base.travel_range_z = [];

p.stage.camera.vendor = '';
p.stage.camera.type = 'tripod';
p.stage.camera.travel_range_z = [];

p.stage.sample.vendor = 'Space Fab';

%% Output
disp( p.stage.rotation )

fprintf( '\nCoherence properties:' )
fprintf( '\n distance source sample: %g', p.coherence.dist_source_sample )
fprintf( '\n FWHM: [%g %g] micron', p.coherence.FWHM_h * 1e6, p.coherence.FWHM_v * 1e6)

fprintf( '\nOptics:\n' )
disp( p.opt )

fprintf( '\nCCD camera:\n' )
disp( p.camera.ccd )

for nn = 1:numel(p.camera.cmos)
    fprintf( '\nCMOS camera %u:\n', nn )
    disp( p.camera.cmos(nn) )
end

fprintf( '\n' )