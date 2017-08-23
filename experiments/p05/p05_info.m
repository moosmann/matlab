% P05 paramter
clear all

%% Optics
p.opt.magn = [5 10 20 40];
p.opt.producer = 'Jena Optik';

%% CCD camera
p.ccd.name = 'EHD SciCam SC09000';
p.ccd.sensor = 'KAF-09000';
p.ccd.dynamic = '16 bit';
p.ccd.pixels = [3056 3056];
p.ccd.pixel_size = 12e-6;
p.ccd.effective_pixel_size_micron = p.ccd.pixel_size * 1e6 ./ p.opt.magn;
p.ccd.pixel_size_micron = p.ccd.pixel_size * 1e6;
p.ccd.fov = p.ccd.pixel_size * p.ccd.pixels(1) ./ p.opt.magn;
p.ccd.fov_mm = p.ccd.fov * 1e3;

%% CMOS KIT camera
p.cmos.name = 'KIT';
p.cmos.sensor = 'CMOSIS CMV20000';
p.cmos.dynamic = '12 bit';
p.cmos.framerate_at_12bit = 30;
p.cmos.pixels = [5120 3840];
p.cmos.pixel_size = 6.4e-6;
p.cmos.effective_pixel_size_micron = p.cmos.pixel_size * 1e6 ./ p.opt.magn;
p.cmos.chip_size_mm = p.cmos.pixels .* p.cmos.pixel_size * 1e3;
p.cmos.pixel_size_micron = p.cmos.pixel_size * 1e6;
p.cmos.fov = p.cmos.pixel_size * p.cmos.pixels(1) ./ p.opt.magn;
p.cmos.fov_mm = p.cmos.fov * 1e3;

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
fprintf( '\nOptics:\n' )
disp( p.opt )

fprintf( '\nCCD camera:\n' )
disp( p.ccd )

fprintf( '\nCMOS camera:\n' )
disp( p.cmos )

disp( p.stage.rotation )