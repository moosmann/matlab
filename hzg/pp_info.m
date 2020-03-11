function p = pp_info( verbose )
% P05 and P07 paramters
if nargin < 1
    verbose = 0;
end

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
p.camera(1).name = 'EHD CCD SciCam SC09000';
p.camera(1).sensor = 'KAF-09000';
p.camera(1).dynamic_bit = 16;
p.camera(1).framerate_16bit_Hz = 0.3;
p.camera(1).pixels = [3056 3056];
p.camera(1).pixel_size_micron = 12;
p.camera(1).effective_pixel_size_micron = p.camera(1).pixel_size_micron ./ p.opt.magn;
p.camera(1).chip_size_mm = p.camera(1).pixels .* p.camera(1).pixel_size_micron * 1e-3;
p.camera(1).fov_5x_mm = 1/5 * p.camera(1).chip_size_mm;
p.camera(1).fov_10x_mm = 1/10 * p.camera(1).chip_size_mm;
p.camera(1).fov_20x_mm = 1/20 * p.camera(1).chip_size_mm;
p.camera(1).fov_40x_mm = 1/40 * p.camera(1).chip_size_mm;

%% CMOS KIT camera
p.camera(2).name = 'KIT';
p.camera(2).sensor = 'CMOSIS CMV20000';
p.camera(2).dynamic_bit = 12;
p.camera(2).framerate_12bit_Hz = 30;
p.camera(2).pixels = [5120 3840];
p.camera(2).pixel_size_micron = 6.4;
p.camera(2).effective_pixel_size_micron = p.camera(2).pixel_size_micron ./ p.opt.magn;
p.camera(2).chip_size_mm = p.camera(2).pixels .* p.camera(2).pixel_size_micron * 1e-3;
p.camera(2).fov_5x_mm = 1/5 * p.camera(2).chip_size_mm;
p.camera(2).fov_10x_mm = 1/10 * p.camera(2).chip_size_mm;
p.camera(2).fov_20x_mm = 1/20 * p.camera(2).chip_size_mm;
p.camera(2).fov_40x_mm = 1/40 * p.camera(2).chip_size_mm;

%% CMOS XIMEA CB500MG
p.camera(3).name = 'XIMEA CB500MG';
p.camera(3).sensor = 'CMOSIS CMV50000';
p.camera(3).dynamic_bit = 12;
p.camera(3).dynamic_range_dB = 64;
p.camera(3).adc_resolution_bit = 12;
p.camera(3).snr_max_dB = 41.6;
p.camera(3).dark_noise_electrons = 8.8;
p.camera(3).dark_current_electron_s = 33;
p.camera(3).dsnu_electrons = 24.5;
p.camera(3).digitization_bit = 12;
p.camera(3).bit_resolutions_bit_pixel = [8,9,10,11,12,16];
p.camera(3).framerate_12bit_Hz = 24;
p.camera(3).framerate_12bit_onchipbin2x_Hz = 30.8;
p.camera(3).expsoure_time_ms = [0.1 1050];
p.camera(3).pixels = [7920 6004];
p.camera(3).pixel_size_micron = 4.6;
p.camera(3).effective_pixel_size_micron = p.camera(3).pixel_size_micron ./ p.opt.magn;
p.camera(3).chip_size_mm = p.camera(3).pixels .* p.camera(3).pixel_size_micron * 1e-3;
p.camera(3).fov_5x_mm = 1/5 * p.camera(3).chip_size_mm;
p.camera(3).fov_10x_mm = 1/10 * p.camera(3).chip_size_mm;
p.camera(3).fov_20x_mm = 1/20 * p.camera(3).chip_size_mm;
p.camera(3).fov_40x_mm = 1/40 * p.camera(3).chip_size_mm;

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
if verbose
    disp( p.stage.rotation )
    
    fprintf( '\nCoherence properties:' )
    fprintf( '\n distance source sample: %g', p.coherence.dist_source_sample )
    fprintf( '\n FWHM: [%g %g] micron', p.coherence.FWHM_h * 1e6, p.coherence.FWHM_v * 1e6)
    
    fprintf( '\nOptics:\n' )
    disp( p.opt )
    
    for nn = 1:numel(p.camera)
        fprintf( '\ncamera %u:\n', nn )
        disp( p.camera(nn) )
    end
    
    fprintf( '\n' )
end