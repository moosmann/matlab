% P05 reconstruction pipeline: preprocessing, filtering, phase retrieval,
% tomographic reconstruction, ...
%
% USAGE:
% Edit parameters in PARAMETERS / SETTINGS section below and run script.
%
% HOW TO RUN THE SCRIPT:
% - Editor windows: press 'F5' when focus is in the Editor window
% - Editor tab: click 'Run' in the toolstrip
% - Command Window: type 'p05_reco' and hit Enter
%
% HOW TO AUTOMATICALLY LOOP RECO OVER DATA SETS:
% To loop over different data or parameters sets see
% 'p05_reco_loop_template' and/or 'p05_create_reco_loop_script'.
%
% For additional information see 'p05_reco_NOTES'.
%
% Please cite following article in the case of publication:
% Moosmann, J. et al. Time-lapse X-ray phase-contrast microtomography for
% in vivo imaging and analysis of morphogenesis Nat. Protocols 9, 294-304
% (2014)
% Also cite the ASTRA Toolbox, see http://www.astra-toolbox.com/
%
% Latest version: https://github.com/moosmann/matlab.git
%
% Written by Julian Moosmann. 

if ~exist( 'external_parameter' ,'var')
    clearvars
end
close all hidden % close all open windows
%dbstop if error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! FAST RECO MODE PARAMTERS !!! OVERWRITES SOME PARAMETERS SET BELOW !!!
fast_reco.run = 1; 
fast_reco.raw_bin = 4;
fast_reco.raw_roi = [0.4 0.6];
fast_reco.proj_range = 4;
fast_reco.ref_range = 10;
%fast_reco.image_correlation.method = 'none';
fast_reco.write.to_scratch = 1;
fast_reco.pixel_filter_radius  = [3 3];
% END OF FAST MODE PARAMTER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scan_path = pwd; % string/pwd. pwd: change to directory of the scan to be reconstructed, sting: absolute scan path
read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
read_flatcor_trafo = @(im) fliplr( im ); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
read_sino_folder = ''; % subfolder to scan path
read_sino_trafo = @(x) rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
pixel_scaling = 1; % to account for beam divergence if pixel size was determined (via MTF) at the wrong distance
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = []; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
% []: use full image;
% [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2); if raw_roi(2) < 0 reads until end - |raw_roi(2)|; relative indexing similar.
% [y0 y1 x0 x1]: vertical + horzontal ROI, each ROI as above
% if -1: auto roi, selects vertical ROI automatically. Use only for DCM. Not working for *.raw data where images are flipped and DMM data.
% if < -1: Threshold is set as min(proj(:,:,[1 end])) + abs(raw_roi)*median(dark(:)). raw_roi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
raw_bin = 2; % projection binning factor: integer
im_trafo = '' ;%'rot90(im,-1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
par.crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
par.stitch_projections = 0; % for 2 pi scans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
par.stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [3 3]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
ring_current_normalization = 0; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
% Correlation of projections and flat fields. Essential for DCM data. Typically improves reconstruction quality of DMM data, too.
% Available methods ('ssim-ml'/'entropy' usually work best):
% 'none' : no correlation, uses median flat, for fast recos
% 'ssim-ml' : Matlab's structural similarity index (SSIM), includes Gaussian smoothing
% 'entropy' : entropy measure of proj over flat, usually similar result to SSIM, but faster
% 'ssim' : own implementation of SSIM without smoothing, usually worse, but sometimes better than 'ssim-ml'
% 'ssim-g' : 'ssim' with smoothing (Gaussian blurring)
% 'cov' : cross covariance
% 'corr' : cross correlation = normalized cross covariance
% 'std' : standard deviation of proj over flat
% 'diff1/2-l1/2': L1/L2-norm of anisotropic (diff1-l*) or isotropic (diff2-l*) difference of projections and flat fields
% 'cross-entropy-*' : asymmetric (12,21) and symmetric (x) cross entropy
image_correlation.force_calc = 0; % bool. force compuation of correlation even though a (previously computed) corrlation matrix exists
image_correlation.num_flats = 3; % number of best maching flat fields used for correction
image_correlation.area_width = [1 100];%[-100 1];% correlation area: index vector or relative/absolute position of [first pix, last pix], negative indexing is supported
image_correlation.area_height = [0.3 0.7]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
ring_filter.apply = 0; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft.dec_levels = 2:5; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft.wname = 'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter.jm.median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 1; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'qp';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 1.7; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, use to determine rot axis positions without processing the full tomogram;
tomo.reco_mode =  '3D';'slice'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = []; %[-.5 .5 -.5 .5 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle.full_range = []; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle.offset = pi; % global rotation of reconstructed volume
tomo.rot_axis.offset = []; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis.position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis.offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.rot_axis.tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
tomo.rot_axis.corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis 
tomo.rot_axis.corr_area2 = []; % ROI to correlate projections at angles 0 & pi 
tomo.fbp_filter.type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
tomo.fbp_filter.freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
tomo.fbp_filter.padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
tomo.fbp_filter.padding_method = 'symmetric';
tomo.butterworth_filter.apply = 0; % use butterworth filter in addition to FBP filter
tomo.butterworth_filter.order = 1;
tomo.butterworth_filter.frequ_cutoff = 0.9;
tomo.astra_pixel_size = 1; % detector pixel size for reconstruction: if different from one 'tomo.vol_size' must to be ajusted, too!
tomo.take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
tomo.algorithm = 'fbp';'sirt'; 'cgls';'sart';'em';'fbp-astra'; % SART/EM only work for 3D reco mode
tomo.iterations = 40; % for 'sirt' or 'cgls'.
tomo.sirt.MinConstraint = []; % If specified, all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions, for example.
tomo.sirt.MaxConstraint = []; % If specified, all values above MaxConstraint will be set to MaxConstraint.
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID, mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder.flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder.phase_map = ''; % subfolder in 'phase_map'
write.subfolder.sino = ''; % subfolder in 'sino'
write.subfolder.reco = ''; % subfolder in 'reco'
write.flatcor = 1; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression.method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression.method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 0; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression.method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression.method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression.method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression.parameter = [0.02 0.02]; % compression-method specific parameter
% dynamic range is compressed s.t. new dynamic range assumes
% 'outlier' : [LOW, HIGH] = write.compression.parameter, eg. [0.01 0.03], outlier given in percent, if scalear LOW = HIGH.
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression.parameter, eg. [-0.01 1]
% 'std' : NUM = write.compression.parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression.parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 0; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 1; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
interactive_mode.angles = 0; % reconstruct slices with different scalings of angles
interactive_mode.angle_scaling_default_search_range = []; % if empty: use a variaton of -/+5 * (angle increment / maximum angle)
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 1; % Interactive retrieval to determine regularization parameter
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.use_cluster = 1; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.poolsize = 0.6; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.7; % Relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
%%% EXPERIMENTAL OR NOT YET IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
verbose = 1;
vert_shift = 0;
offset_shift = 0;
scan_position = [];
logpar = [];
                
%%% Parameters set by reconstruction loop script 'p05_reco_loop' %%%%%%%%%%
if exist( 'external_parameter' ,'var')
    par.visual_output = 0;
    dbstop if error
    interactive_mode.rot_axis_pos = 0;
    fast_reco = 0;
    fn = fieldnames( external_parameter );
    for nn = 1:numel( fn )
        var_name = fn{nn};
        var_val = getfield( external_parameter, var_name );
        assignin('caller', var_name, var_val )
    end
    clear external_parameter;
    dbstop if error
end

%% % FAST RECO MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist( 'fast_reco', 'var') && fast_reco.run
    fn = fieldnames( fast_reco );
    for nn = 1:numel( fn )
        name = fn{nn};
        if ~strcmp( name, 'run' )
           val = fast_reco.(name);
           assign_variable( name, val );
        end
    end
    cprintf( 'Red', '\nATTENTION: fast reco mode is turned on!\n\n' )
end
%%

% Parameter checks
if ~phase_retrieval.apply
    phase_retrieval.post_binning_factor = 0;
end
if interactive_mode.rot_axis_tilt && strcmpi( tomo.reco_mode, 'slice' )
    error( 'Slicewise reconstruction and reconstruction with tilted rotation axis are not compatible!' )
end
% infer (relative) vol_shape from vol_size if vol_shape is empty, but vol_size is given
if ~isempty( tomo.vol_size ) && isempty( tomo.vol_shape )
    tomo.vol_shape = tomo.vol_size(2:2:end) - tomo.vol_size(1:2:end);
end
if ~isempty( tomo.rot_axis.offset ) && ~isempty( tomo.rot_axis.position )
    error('tomo.rot_axis.offset (%f) and tomo.rot_axis.position (%f) cannot be used simultaneously. One must be empty.', tomo.rot_axis.offset, tomo.rot_axis.position)
end

% Default assignment if non-existing or empty!
assign_default( 'tomo.rot_axis.offset', 0 )
assign_default( 'pixel_filter_radius', [3 3] )
assign_default( 'image_correlation.force_calc', 0 );
assign_default( 'write.path', '' )
assign_default( 'write.parfolder', '' )
assign_default( 'write.subfolder.reco', '' )
assign_default( 'write.subfolder.flatcor', '' )
assign_default( 'write.subfolder.phase_map', '' )
assign_default( 'write.subfolder.sino', '' )
assign_default( 'write.sino_shift_cropped', 0 ) 
assign_default( 'write.deleteFiles', 0)
assign_default( 'write.beamtimeID', '' )
assign_default( 'tomo.reco_mode', '3D' )
assign_default( 'tomo.rot_axis.tilt', 0 )
assign_default( 'write.scan_name_appendix', '' )
assign_default( 'interactive_mode.rot_axis_pos_default_search_range', -4:0.5:4 ) % binned pixel
assign_default( 'interactive_mode.rot_axis_tilt_default_search_range', -0.005:0.001:0.005 ) % radian
assign_default( 'interactive_mode.phase_retrieval_default_search_range', [] )
assign_default( 'interactive_mode.angles', 0 );
assign_default( 'interactive_mode.angle_scaling_default_search_range', [] );

assign_default( 'tomo.rot_axis.corr_area2', [0.1 0.9] );
 
% Define variables from struct fields for convenience
raw_bin = single( raw_bin );
par.raw_bin = raw_bin;
phase_bin = phase_retrieval.post_binning_factor; % alias for readablity
eff_pixel_size_binned = raw_bin * eff_pixel_size;

astra_clear % if reco was aborted, ASTRA memory is not cleared

% Utility functions
imsc1 = @(im) imsc( rot90( im ) );
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';

% Disable warnings
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag' );
warning( 'off', 'MATLAB:hg:AutoSoftwareOpenGL' );

fprintf( 'START RECONSTRUCTION: ')

%% Folders

% Scan path
while scan_path(end) == filesep
    scan_path(end) = [];
end
[raw_path, scan_name] = fileparts(scan_path);
scan_path = [scan_path, filesep];
[beamtime_path, raw_folder] = fileparts(raw_path);
[~, beamtime_id] = fileparts(beamtime_path);
if ~strcmp(raw_folder, 'raw') && ~read_sino && ~read_flatcor
    error('Given path does not contain a ''raw'' folder: %s', raw_folder)
end

% Output path and folder
out_folder = 'processed';
if write.to_scratch
    out_folder = 'scratch_cc';
end
if isempty( write.path )
    write.path = [beamtime_path, filesep, out_folder, filesep, scan_name];
else
    write.path = [write.path, filesep, scan_name];
end
if ~isempty( write.scan_name_appendix )
    write.path = [ write.path '_' write.scan_name_appendix ];
end
write.parpath = [write.path filesep ];
if ~isempty(write.parfolder)
    write.path = [write.path, filesep, write.parfolder];
end

% Save raw path to file for shell short cut
filename = [userpath, filesep, 'experiments/p05/path_to_latest_raw'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', raw_path );
fclose( fid );

fprintf( '%s', scan_name )
fprintf( ' at %s', datetime )
fprintf( '\n scan_path:\n  %s', scan_path )

% Memory
fprintf( '\n user :  %s', getenv( 'USER' ) );
fprintf( '\n hostname : %s', getenv( 'HOSTNAME' ) );
[mem_free, mem_avail, mem_total] = free_memory;
fprintf( '\n system memory: free, available, total : %.0f GiB, %.0f GiB, %.0f GiB', round([mem_free/1024^3, mem_avail/1024^3, mem_total/1024^3]) )
if isempty( tomo.astra_gpu_index )
    tomo.astra_gpu_index = 1:gpuDeviceCount;
end

% GPU info
fprintf( '\n GPU : [index, total memory/GiB] =' )
for mm = 1:numel( tomo.astra_gpu_index )
    nn = tomo.astra_gpu_index(mm);
    gpu = parallel.gpu.GPUDevice.getDevice( nn );
    mem_total = gpu.TotalMemory/1024^3;
    fprintf( ' [%u %.3g]', nn, mem_total )
end

% Save scan path to file
filename = [userpath, filesep, 'experiments/p05/path_to_latest_scan'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', scan_path );
fclose( fid );

% Reco path
if isempty( write.subfolder.reco )
    reco_path = [write.path, filesep, 'reco', filesep];
else
    reco_path = [write.path, filesep, 'reco', filesep, write.subfolder.reco, filesep];
end
write.reco_path = reco_path;

% Figure path
write.fig_path = [write.path, filesep, 'figures', filesep];
fig_path = write.fig_path;
CheckAndMakePath( fig_path )

if ~read_flatcor && ~read_sino
    
    % Wait until scan finishes
    filename = sprintf( '%s%sscan.log', scan_path, scan_name );
    if ~exist( filename, 'file' )
        % back up for older log file name schemes
        str = dir( sprintf( '%s*scan.log', scan_path) );
        filename = sprintf( '%s/%s', str.folder, str.name);
    end
    while exist(filename, 'file' ) && ~getfield( dir( filename ) ,'bytes')
        fprintf( '\nWaiting for scan to finish.' )
        pause(1);
    end
    
    % Path to flat-field corrected projections
    flatcor_path = sprintf( '%s/flat_corrected/rawBin%u/', write.path, raw_bin );
    if ~isempty( write.subfolder.flatcor )
        flatcor_path =  sprintf( '%s%s/', flatcor_path, write.subfolder.flatcor );
    end
    PrintVerbose( write.flatcor, '\n flatcor_path:\n  %s', flatcor_path)
    
    % Path to retrieved phase maps
    write.phase_map_path = sprintf( '%s/phase_map/rawBin%u/', write.path, raw_bin );
    if ~isempty( write.subfolder.phase_map )
        %write.phase_map_path = [write.path, filesep, 'phase_map', filesep, write.subfolder.phase_map, filesep];
        write.phase_map_path = sprintf( '%s%s/',  write.phase_map_path, write.subfolder.phase_map );
    end
    PrintVerbose( write.phase_map, '\n phase_map_path:\n  %s', write.phase_map_path)
    
    % Sinogram path
    sino_path = sprintf( '%s/sino/rawBin%u/', write.path, raw_bin );
    write.sino_phase_path = sprintf( '%s/sino_phase/rawBin%u/', write.path, raw_bin );
    if ~isempty( write.subfolder.sino )
        sino_path = sprintf( '%s%s/', sino_path, write.subfolder.sino );
        write.sino_phase_path = sprintf( '%s%s/', write.sino_phase_path, write.subfolder.sino );
    end
    PrintVerbose( write.sino, '\n sino_path:\n  %s', sino_path)
    PrintVerbose( phase_retrieval.apply & write.phase_sino, '\n sino_phase_path:\n  %s', write.sino_phase_path)
    
    % Projection file names
    proj_names = FilenameCell( [scan_path, '*.img'] );
    raw_data = 0;
    if isempty( proj_names )
        proj_names =  FilenameCell( [scan_path, '*img*.tif'] );
        raw_data = 0;
    end
    if isempty( proj_names )
        proj_names =  FilenameCell( [scan_path, '*img*.raw'] );
        raw_data = 1;
    end
    if isempty( proj_names )
        proj_names =  FilenameCell( [scan_path, '*proj*.tif'] );
        raw_data = 0;
    end
    if isempty( proj_names )
        proj_names =  FilenameCell( [scan_path, '*proj*.raw'] );
        raw_data = 1;
    end
    num_proj_found = numel(proj_names);
    
    % Ref file names
    ref_names = FilenameCell( [scan_path, '*.ref'] );
    if isempty( ref_names )
        ref_names = FilenameCell( [scan_path, '*ref.tif'] ); 
    end
    if isempty( ref_names )
        ref_names =  FilenameCell( [scan_path, '*flat*.tif'] );
    end
    if isempty( ref_names )
        ref_names =  FilenameCell( [scan_path, '*ref*.raw'] );
    end
    if isempty( ref_names )
        ref_names =  FilenameCell( [scan_path, '*flat*.raw'] );
    end
    if isempty( ref_names )
        ref_names = FilenameCell( [scan_path, '*ref*.tif'] ); 
    end
    num_ref_found = numel(ref_names);
    if isempty( ref_range )
        ref_range = 1;
    end
    if numel( ref_range ) == 1
        ref_range = 1:ref_range:num_ref_found;
    end
    % position of running index
    re = regexp( ref_names{1}, '\d{6,6}');
    if numel( re ) == 1
        imtype_str_flag = re;
    elseif strcmpi( ref_names{1}(end-6:end-4), 'ref' )
        imtype_str_flag = '64'; % 0
    elseif strcmpi( ref_names{1}(end-11:end-9), 'ref' )
        imtype_str_flag = '119'; % 1
    end
    num_ref_used = numel( ref_range );
    ref_names_mat = NameCellToMat( ref_names(ref_range) );
    ref_nums = CellString2Vec( ref_names(ref_range) , imtype_str_flag);
    fprintf( '\n refs found : %g', num_ref_found)
    fprintf( '\n refs used : %g', num_ref_used)
    fprintf( '\n reference range used : %g:%g:%g%', ref_range(1), ref_range(2) - ref_range(1), ref_range(end))
    
    % Dark file names
    dark_names = FilenameCell( [scan_path, '*.dar'] );
    if isempty( dark_names )
        dark_names = FilenameCell( [scan_path, '*dar.tif'] );
    end    
    if isempty( dark_names )
        dark_names =  FilenameCell( [scan_path, '*dar*.raw'] );
    end
    if isempty( dark_names )
        dark_names = FilenameCell( [scan_path, '*dar*.tif'] );
    end
    
    dark_nums = CellString2Vec( dark_names, imtype_str_flag );
    num_dark = numel(dark_names);
    fprintf( '\n darks found : %g', num_dark)
    
    % Projection range to read
    if isempty( proj_range )
        proj_range = 1;
    end
    if numel( proj_range ) == 1
        proj_range = 1:proj_range:num_proj_found;
    end
    num_proj_used = numel( proj_range );
    proj_nums = CellString2Vec( proj_names(proj_range), imtype_str_flag);
    fprintf( '\n projections found : %g', num_proj_found)
    fprintf( '\n projections used : %g', num_proj_used)
    fprintf( '\n projection range used : first:stride:last =  %g:%g:%g', proj_range(1), proj_range(2) - proj_range(1), proj_range(end))
    
    %% Log files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hdf5 log from statussever
    h5log = dir( sprintf('%s*_nexus.h5', scan_path) );
    if numel( h5log ) == 1
        h5log = [h5log.folder filesep h5log.name];
    end
    im_shape_raw = [];
    dtype = '';
    tif_info = [];
    % old scan-log, still needed
    filename = sprintf( '%s%sscan.log', scan_path, scan_name );
    if ~exist( filename, 'file' )
        % back up for older log file name schemes
        str = dir( sprintf( '%s*scan.log', scan_path) );
        filename = sprintf( '%s/%s', str.folder, str.name);
    end
    [logpar, cur, cam] = p05_log( filename );
    if isempty( eff_pixel_size )
        eff_pixel_size = logpar.eff_pixel_size;
    end
    % Scaling of pixel size if MTF is wrong. Important for lateral shift scans
    if exist( 'pixel_scaling', 'var' ) && ~isempty( pixel_scaling )
        eff_pixel_size = pixel_scaling * eff_pixel_size;
    end
    eff_pixel_size_binned = raw_bin * eff_pixel_size;
    exposure_time = logpar.exposure_time;
    if isempty( energy )
        if isfield( logpar, 'energy')
            energy = logpar.energy;
        end
    end
    if isempty( sample_detector_distance )
        sample_detector_distance = logpar.sample_detector_distance;
    end
    if ~exist( h5log, 'file')
        % Image shape and ROI
        filename = sprintf('%s%s', scan_path, ref_names{1});
        if ~raw_data
            [im_raw, tif_info] = read_image( filename, '', [], [], [], '', im_trafo );
            im_shape_raw = size( im_raw );
        else
            switch lower( cam )
                case 'ehd'
                    im_shape_raw = [3056 3056];
                    dtype = 'uint16';
                case 'kit'
                    im_shape_raw = [5120 3840];
                    dtype = 'uint16';
            end
            eff_pixel_size_binned = raw_bin * eff_pixel_size;
            im_raw = read_raw( filename, im_shape_raw, dtype );
        end
    else
        % HDF5 log
        h5log_info = h5info( h5log );
        % energy, exposure time, image shape
        switch lower( cam )
            %%% CHECK h5 entry of camera1 / camera2 !!
            case 'ehd'
                if isempty( energy )
                    energy = double(h5read( h5log, '/entry/hardware/camera1/calibration/energy') );
                end
                exposure_time = double(h5read( h5log, '/entry/hardware/camera1/calibration/exptime') );
                %im_shape_raw = [3056 3056];
                dtype = 'uint16';
            case 'kit'
                if isempty( energy )
                    energy = double(h5read( h5log, '/entry/hardware/camera2/calibration/energy') );
                end
                exposure_time = double(h5read( h5log, '/entry/hardware/camera2/calibration/exptime') );
                %im_shape_raw = [5120 3840];
                dtype = 'uint16';
        end
        
        eff_pixel_size_binned = raw_bin * eff_pixel_size;
        % Image shape
        filename = sprintf('%s%s', scan_path, ref_names{1});
        
        % mod: breaks raw data support
        % Fixed: 2019-07-10
        % CLEAN UP required
        if raw_data
            switch lower( cam )
                case 'ehd'
                    im_shape_raw = [3056 3056];
                    %dtype = 'uint16';
                case 'kit'
                    im_shape_raw = [5120 3840];
                    %dtype = 'uint16';
            end
            [im_raw, tif_info] = read_image( filename, '', [], tif_info, im_shape_raw, dtype, im_trafo );
        else
            [im_raw, tif_info] = read_image( filename, '', [], tif_info, [], dtype, im_trafo );
        end
        im_shape_raw = size( im_raw );
        
        % images
        stimg_name.value = unique( h5read( h5log, '/entry/scan/data/image_file/value') );
        stimg_name.time = h5read( h5log,'/entry/scan/data/image_file/time');
        stimg_key.value = h5read( h5log,'/entry/scan/data/image_key/value');
        stimg_key.time = double( h5read( h5log,'/entry/scan/data/image_key/time') );
        % PETRA ring current
        [petra.time, index] = unique( h5read( h5log,'/entry/hardware/beam_current/current/time') );
        petra.current = h5read( h5log,'/entry/hardware/beam_current/current/value');
        petra.current = petra.current(index);
        if par.visual_output
            name = 'PETRA beam current from status server';
            f = figure( 'Name', name, 'WindowState', 'maximized');
            x = double(petra.time(2:1:end)-petra.time(2)) / 1000 / 60;
            y = petra.current(2:1:end);
            plot( x, y, '.' )
            xlabel( 'time / min' )
            ylabel( 'current / mA' )
            title( name )
            axis tight
            drawnow
            CheckAndMakePath( fig_path )
            fig_filename = sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) ;
            saveas( f, fig_filename );
        end
        
        % rotation axis
        s_rot.time = double( h5read( h5log, '/entry/scan/data/s_rot/time') );
        s_rot.value = h5read( h5log, '/entry/scan/data/s_rot/value');
        
        %% Lateral shift
        h5log_group = h5info(h5log, '/entry/scan/data/' );
        if sum( strcmp('/entry/scan/data/s_stage_x',{h5log_group.Groups.Name}))
            % Read out lateral rotation axis shift form log file
            s_stage_x.time = double( h5read( h5log, '/entry/scan/data/s_stage_x/time') );
            s_stage_x.value = h5read( h5log, '/entry/scan/data/s_stage_x/value');
            
            if numel( s_stage_x.value )
                offset_shift_micron = s_stage_x.value( ~boolean( stimg_key.value(logpar.n_dark+1:end) ) );
                % Shift or static position
                if std( offset_shift_micron )
                    offset_shift_micron = offset_shift_micron(proj_range);
                    offset_shift = 1e-3 / eff_pixel_size * offset_shift_micron;
                    offset_shift = 1 + offset_shift - min( offset_shift(:) ) ;
                    
                    % Overwrite lateral shift if offset shift is provided as parameter
                    if isequal( std( offset_shift ), 0 ) && ~isempty( tomo.rot_axis.offset_shift )
                        offset_shift = tomo.rot_axis.offset_shift / raw_bin * (0:num_proj_used) / num_proj_used;
                    end
                    
                    % Tranform to integer pixel-wise shifts
                    tmp = offset_shift;
                    offset_shift = round( offset_shift );
                    if std( offset_shift - tmp ) > 1e-2
                        error( 'Offset shift not on integer pixel scale' )
                    end
                    
                    % Plot offset shift
                    if par.visual_output
                        f = figure( 'Name', 'rotation axis offset shift', 'WindowState', 'maximized');
                        plot( offset_shift, '.')
                        title( sprintf('Rotation axis offset shift') )
                        axis equal tight
                        xlabel( 'projection number' )
                        ylabel( 'lateral shift / pixel' )
                        legend( sprintf( 'effective pixel size binned: %.2f micron', eff_pixel_size_binned * 1e6 ) )
                        drawnow
                        CheckAndMakePath( fig_path )
                        saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
                    end
                    
                    % Lateral scanning
                    % Position index extracted by jump in offset_shift
                    scan_position_index = zeros( size( offset_shift ) );
                    pos = 1;
                    scan_position_index(1) = pos;
                    for nn = 2:numel( offset_shift )
                        if abs( offset_shift(nn) - offset_shift(nn-1) ) > 201
                            pos = pos + 1;
                        end
                        scan_position_index(nn) = pos;
                    end
                    % Absolute scan position without lateral offset
                    scan_position = zeros( size( offset_shift ) );
                    for nn = 1:pos
                        m = scan_position_index == nn;
                        scan_position(m) = min( offset_shift(m) ) - 1;
                    end
                    offset_shift = offset_shift - scan_position;
                    
                    % Scale position because of binning for tomo reco
                    scan_position = scan_position + mean( scan_position );
                    scan_position = 1 / raw_bin * scan_position;
                    
%                     Y = [ normat( scan_position_index ), normat(offset_shift )];
%                     plot( Y, '.' )
%                     axis tight

                end % if std( offset_shift_micron )
            end % if numel( s_stage_x.value )
        end
        
        %% Vertical shift
        if sum( strcmp('/entry/scan/data/s_stage_z',{h5log_group.Groups.Name}))
            % Read out vertical shift form HDF5 log file
            s_stage_z.time = double( h5read( h5log, '/entry/scan/data/s_stage_z/time') );
            s_stage_z.value = h5read( h5log, '/entry/scan/data/s_stage_z/value');
            
            % Static or shift?
            if numel( s_stage_z.value ) % is not empty
                switch numel( s_stage_z.value )
                    case num_proj_found
                         vert_shift_micron = s_stage_z.value;
                    case num_proj_found + num_dark + num_ref_found
                        m = stimg_key.value == 0 ;
                        vert_shift_micron = s_stage_z.value( m );
                    case num_proj_found + num_ref_found
                        %m = stimg_key.value(logpar.n_dark + 1:end) == 0 ;
                        %m = stimg_key.value(num_proj_found + 1:end) == 0 ;
                        vert_shift_micron = s_stage_z.value( round(num_ref_found/2) + (1:num_proj_found) );
                    otherwise
                        vert_shift_micron = s_stage_z.value( 1:num_proj_found );    
                end
                if std( vert_shift_micron )
                    vert_shift_micron = vert_shift_micron(proj_range);
                    
                    % Check
                    vert_shift = vert_shift_micron * 1e-3 / eff_pixel_size_binned;
                    vert_shift = SubtractMean( vert_shift );                                       
                   
                    fprintf( '\n vertical shift absolute / micron : [%g %g]', min( vert_shift_micron ), max( vert_shift_micron ) )
                    fprintf( '\n vertical shift relative / binned pixel : [%g %g] ', min( vert_shift), max( vert_shift) )
                    dz_micron = max( vert_shift_micron ) - min(vert_shift_micron);
                    dz = max( vert_shift ) - min(vert_shift );
                    fprintf( '\n vertical shift diff : %g micron, %g binned pixel', dz_micron, dz )
                    fprintf( '\n vertical shift / #proj / unbinned pixel : %g', dz / num_proj_found * raw_bin );
                    
                    % Plot vertical shift
                    if par.visual_output
                        name = 'spiral scan: vertical shift';
                        f = figure( 'Name', name, 'WindowState', 'maximized');
                        plot( vert_shift, '.')
                        title( name )
                        axis equal tight
                        xlabel( 'projection number' )
                        ylabel( 'vertical / pixel' )
                        legend( sprintf( 'shift / #proj / pixel %f', dz / num_proj_found ) )
                        drawnow
                        CheckAndMakePath( fig_path )
                        saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
                    end
                end % if std( vert_shift_micron )
            end % if numel( s_stage_z.value )
        end %if sum( strcmp('/entry/scan/data/s_stage_z',{a.Groups.Name}))
        
        %% Ring current
        X = double( petra.time(2:end) ); % first value is zero
        V = double( petra.current(2:end) ); % first value is zero
        Xq = double( stimg_name.time );
        stimg_name.current = (interp1( X, V, Xq, 'next', 100) + interp1( X, V, Xq + exposure_time, 'previous', 100) ) / 2;
        cur_ref_val = stimg_name.current( stimg_key.value == 1 );
        cur_ref_name = stimg_name.value( stimg_key.value == 1 );
        re = regexp( cur_ref_name{1}, '\d{6,6}');
        if numel( re ) == 1
            imtype_str_flag = re;
        elseif strcmpi( ref_names{1}(end-6:end-4), 'ref' )
            imtype_str_flag = '64'; % 0
        elseif strcmpi( ref_names{1}(end-11:end-9), 'ref' )
            imtype_str_flag = '119'; % 1
        end
        for nn = numel( cur_ref_name ):-1:1
            cur.ref(nn).val = cur_ref_val(nn);
            cur.ref(nn).name = cur_ref_name{nn};
            if imtype_str_flag >= 0
                cur.ref(nn).ind = str2double(cur.ref(nn).name(imtype_str_flag + (0:5)));
            elseif imtype_str_flag == -1
                cur.ref(nn).ind = str2double(cur.ref(nn).name(end-12:end-8));
            elseif imtype_str_flag == -2
                cur.ref(nn).ind = str2double(cur.ref(nn).name(end-7:end-4));
            end
        end
        cur_proj_val = stimg_name.current( stimg_key.value == 0);
        cur_proj_name = stimg_name.value( stimg_key.value == 0);
        for nn = numel( cur_proj_name ):-1:1
            cur.proj(nn).val = cur_proj_val(nn);
            cur.proj(nn).name = cur_proj_name{nn};
            if imtype_str_flag >= 0
                cur.proj(nn).ind = str2double(cur.proj(nn).name(imtype_str_flag + (0:5)));
            elseif imtype_str_flag == -1
                cur.proj(nn).ind = str2double(cur.proj(nn).name(end-12:end-8));
            elseif imtype_str_flag == -2
                cur.proj(nn).ind = str2double(cur.proj(nn).name(end-7:end-4));
            end
        end
    end % if ~exist( h5log, 'file')
    
    %% Raw ROI
    raw_roi = set_raw_roi( raw_roi, par, im_shape_raw, im_raw, tif_info, dtype, im_trafo, scan_path, fig_path, ref_names, dark_names );
    
    %% Print info
    im_roi = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo );
    im_shape_roi = size( im_roi );
    im_shape_binned1 = floor( size( im_roi, 1 ) / raw_bin );
    im_shape_binned2 = floor( size( im_roi, 2 ) / raw_bin );
    fprintf( '\n energy : %.1f keV', energy / 1e3 )
    fprintf( '\n distance sample dector : %.1f mm', sample_detector_distance * 1000 )
    fprintf( '\n effective pixel size unbinned : %.2f micron',  eff_pixel_size * 1e6)
    fprintf( '\n effective pixel size binned: %.2f micron',  eff_pixel_size_binned * 1e6)
    fprintf( '\n image size : %.2f mm x %.2f mm', im_shape_raw * eff_pixel_size *1e3 )
    fprintf( '\n image shape : %u x %u = %.1g pixels', im_shape_raw, prod( im_shape_raw ))
    fprintf( '\n image shape roi : %u x %u = %.1g pixels', im_shape_roi, numel( im_roi ) )
    numel_im_roi_binned = im_shape_binned1 * im_shape_binned2;
    fprintf( '\n image shape roi binned : %u x %u = %.1g pixels', im_shape_binned1, im_shape_binned2, numel_im_roi_binned)
    fprintf( '\n raw binning factor : %u', raw_bin)
    
    % GPU info
    for mm = numel( tomo.astra_gpu_index ):-1:1
        nn = tomo.astra_gpu_index(mm);
        gpu(nn) = gpuDevice(nn);
        mem_avail(nn) =  gpu(nn).AvailableMemory;
        mem_total(nn) = gpu(nn).TotalMemory;
        fprintf( '\n GPU %u memory available : %.3g GiB (%.2f%%) of %.3g GiB', nn, mem_avail(nn)/1024^3, 100*mem_avail(nn)/mem_total(nn), mem_total(nn)/1024^3 )
    end
    gpu_mem_requ_per_im = prod( im_shape_roi ) * (4 + 2 + 2 + 1 ); % Limit by GPU pixel filter: 1 single + 2 uint16 + 1 logical
    poolsize_max = floor( par.poolsize_gpu_limit_factor * min( mem_avail ) / gpu_mem_requ_per_im / numel( tomo.astra_gpu_index ) );
    fprintf( ' \n estimated GPU memory required per image for pixel filtering : %g MiB', gpu_mem_requ_per_im / 1024^2 )
    fprintf( ' \n GPU poolsize limit factor : %g', par.poolsize_gpu_limit_factor )
    fprintf( ' \n GPU memory induced maximum poolsize : %u', poolsize_max )
    
    % Start parallel CPU pool
    parpool_tmp_folder = [beamtime_path filesep 'scratch_cc'];
    [poolobj, par.poolsize] = OpenParpool( par.poolsize, par.use_cluster, parpool_tmp_folder, 0, poolsize_max );
    
    %% Dark field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = toc;
    fprintf( '\nProcessing %u dark fields.', num_dark)
    dark = zeros( [im_shape_roi, num_dark], 'uint16');
    fprintf( ' Allocated memory: %.2f MiB,', Bytes( dark, 2 ) )
    filt_pix_par.threshold_hot = pixel_filter_threshold_dark(1);
    filt_pix_par.threshold_dark = pixel_filter_threshold_dark(2);
    filt_pix_par.medfilt_neighboorhood = pixel_filter_radius;
    filt_pix_par.filter_dead_pixel = 1;
    filt_pix_par.filter_Inf = 0;
    filt_pix_par.filter_NaN = 0;
    filt_pix_par.verbose = 0;
    
    %startS = ticBytes( gcp );
    parfor nn = 1:num_dark
        
        % Read image
        filename = sprintf('%s%s', scan_path, dark_names{nn});
        im_int = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo);
        
        % Remove large outliers. Assume Poisson distribtion at large lambda
        % is approximately a Gaussian distribution and set all value above
        % mean + 4 * std (99.994 of values lie within 4 std). Due to
        % outliers 4*std will contain much more values and is a good
        % estimate
        im_float = single( im_int(:) );
        im_mean = mean( im_float );
        im_std = std( im_float );
        im_int( im_int > im_mean + 4*im_std) = uint16( im_mean);
        
        % Filter pixels
        im_int = FilterPixelGPU( im_int, filt_pix_par);
        
        % Assign image to stack
        dark(:, :, nn) = im_int;
    end
    %toc_bytes.read_dark = tocBytes( gcp, startS );
    
    darks_min = min( dark(:) );
    darks_max = max( dark(:) );
    
    % Reject dark images which are all zero
    darks_to_use = zeros( 1, num_dark, 'logical' );
    parfor nn = 1:num_dark
        darks_to_use(nn) = boolean( max2( dark(:,:,nn) )  );
    end
    
    % Median dark
    dark = squeeze( median(dark(:,:,darks_to_use), 3) );
    dark_med_min = min( dark(:) );
    dark_med_max = max( dark(:) );
    
    % Binned dark
    dark_binned = 1 / raw_bin^2 * Binning( dark, raw_bin);
    fprintf( ' done in %.1f s', toc-t)
    fprintf( '\n min/max of all darks : %g %g', darks_min, darks_max);
    fprintf( '\n min/max of median dark : %g %g', dark_med_min, dark_med_max);
    
    % Fig: raw + dark field
    if par.visual_output
        h1 = figure( 'Name', 'data and flat-and-dark-field correction', 'WindowState', 'maximized');
        subplot(2,3,1)
        imsc1( dark_binned );
        title(sprintf('median dark field'))
        colorbar
        axis equal tight
        xticks('auto'),yticks('auto')
        drawnow
    end
    
    %% Flat field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = toc;
    fprintf( '\nProcessing %u flat fields.', num_ref_used )
    % Correlation roi area
    flat_corr_area1 = image_correlation.area_width;
    flat_corr_area2 = image_correlation.area_height;
    flat_corr_area1 = IndexParameterToRange( flat_corr_area1, im_shape_roi(1) );
    flat_corr_area2 = IndexParameterToRange( flat_corr_area2, im_shape_roi(2) );
    % Preallocation
    flat = zeros( [im_shape_binned1, im_shape_binned2, num_ref_used], 'single');
    roi_flat = zeros( numel( flat_corr_area1 ), numel( flat_corr_area2 ) , num_proj_used, 'single');
    num_zeros = zeros( 1, num_ref_used );
    fprintf( ' Allocated memory: %.2f GiB,', Bytes( flat, 3 ) )
    refs_to_use = zeros( 1, size( flat,3), 'logical');
    filt_pix_par.threshold_hot = pixel_filter_threshold_flat(1);
    filt_pix_par.threshold_dark = pixel_filter_threshold_flat(2);
    
    % Parallel loop over refs
    startS = ticBytes( gcp );
    parfor nn = 1:num_ref_used
               
        % Read
        filename = sprintf('%s%s', scan_path, ref_names_mat(nn,:));
        im_int = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo );
        
        % Filter pixel
        im_int = FilterPixelGPU( im_int, filt_pix_par);
        
        % Dark field correction
        im_int = im_int - dark;
        
        % Correlation ROI
        roi_flat(:,:,nn) = im_int(flat_corr_area1,flat_corr_area2);
        
        % Binning
        im_float_binned = Binning( im_int, raw_bin) / raw_bin^2;        
        
        % Count for zeros
        num_zeros(nn) =  sum( im_float_binned(:) < 1  );
        
        % Discard if any pixel is zero.
        refs_to_use(nn) = ~boolean( num_zeros(nn)  );
        
        % Assign image to stack
        flat(:, :, nn) = im_float_binned;
    end
    toc_bytes.read_flat = tocBytes( gcp, startS );
    
    % Delete empty refs
    zz = ~refs_to_use;
    if sum( zz(:) )
        fprintf( ' Deleting empty flats' )
        flat(:,:,~refs_to_use) = [];
    end
    
    % min/max values before dark field subtraction and ring current normalization
    flat_min = min( flat(:) );
    flat_max = max( flat(:) );
    
    % Ring current normalization
    if ring_current_normalization(1)
        ref_ind_from_filenames = ref_nums;
        ref_ind_from_log = [cur.ref(ref_range).ind];
        if isequal( ref_ind_from_filenames, ref_ind_from_log )
            ref_rc = [cur.ref(ref_range).val];
            ref_rcm = mean( ref_rc(:) );
            scale_factor = 100 ./ shiftdim( ref_rc(refs_to_use), -1 );
            flat = bsxfun( @times, flat, scale_factor );
            if par.visual_output
                hrc = figure( 'Name', 'PETRA III beam current: Interpolation at image time stamps', 'WindowState', 'maximized');
                subplot(1,1,1);
                plot( ref_rc(:), '.' )
                axis tight
                title(sprintf('ring current: flat fields'))
                legend( sprintf( 'mean: %.2f mA', ref_rcm) )
                drawnow
            end
        else
            fprintf('\n WARNING: flat fields not normalized by ring current. Names read from directory and log-file are inconsistent.\n')
        end
    end
    
    flat_min2 = min( flat(:) );
    flat_max2 = max( flat(:) );
    nn =  sum( flat(:) < 1 );
    if nn > 0
        fprintf('\n WARNING: flat field contains %u zeros\n', nn)
    end
    
    nn = sum( ~refs_to_use(:) );
    num_ref_used = num_ref_used - nn;
    fprintf( ' done in %.1f s', toc-t)
    fprintf( '\n min/max of all flats : %6g %6g', flat_min, flat_max);
    fprintf( '\n min/max of all corrected flats : %6g %6g', flat_min2, flat_max2);
    
    PrintVerbose( nn,'\n discarded empty refs : %u, %.2f%%', nn, 100 * nn / num_ref_found )
    if sum( num_zeros )
        fprintf( '\n flat fields with zeros :' )
        % print #zeros if not all pixels are zero
        for nn = 1:numel(num_zeros)
            if num_zeros(nn) ~= 0
                if isequal( num_zeros(nn), numel_im_roi_binned )
                    fprintf( ' %u', nn )
                else
                    fprintf( ' %u:%u', nn, num_zeros(nn) )
                end
            end
        end
    end
    
    % Show flat field
    if par.visual_output
        if exist( 'h1' , 'var' ) && isvalid( h1 )
            figure(h1)
        else
            h1 = figure( 'Name', 'data and flat-and-dark-field correction', 'WindowState', 'maximized');
        end
        subplot(2,3,2)
        imsc1( flat(:,:,1) )
        title(sprintf('flat field #1'))
        colorbar
        axis equal tight
        drawnow
    end
    
    %% Projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = toc;
    fprintf( '\nProcessing %u projections.', num_proj_used )
    img_names_mat = NameCellToMat( proj_names(proj_range) );
    
    filt_pix_par.threshold_hot = pixel_filter_threshold_proj(1);
    filt_pix_par.threshold_dark = pixel_filter_threshold_proj(2);
    
    % Display first raw image
    if par.visual_output
        if exist( 'h1' , 'var' ) && isvalid( h1 )
            figure(h1)
        else
            h1 = figure( 'Name', 'data and flat-and-dark-field correction', 'WindowState', 'maximized');
        end
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        im_int = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo );
        %im_int = FilterPixel( im_int, pixel_filter_threshold_proj, 0, pixel_filter_radius);
        im_int = FilterPixelGPU( im_int, filt_pix_par );
        raw1 = Binning( im_int, raw_bin) / raw_bin^2;
        subplot(2,3,3)
        imsc1( raw1 )
        title(sprintf('raw proj #1'))
        colorbar
        axis equal tight
        drawnow
    end
    
    % Get absolut filter thresholds from percentage-wise pixel filtering
    % of 1st, middle, and last projection to speed up processing
    if filt_pix_par.threshold_hot  < 1 || filt_pix_par.threshold_dark < 0.5
        filename = sprintf('%s%s', scan_path, img_names_mat(num_proj_used, :));
        im_int = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo );
        [~, ht(3), dt(3)] = FilterPixelGPU( im_int, filt_pix_par );
        %[~, ht(3), dt(3)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo ), pixel_filter_threshold_proj, 0, pixel_filter_radius);
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        im_int = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo );
        [~, ht(2), dt(2)] = FilterPixelGPU( im_int, filt_pix_par );
        %[~, ht(2), dt(2)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo ), pixel_filter_threshold_proj, 0, pixel_filter_radius);
        filename = sprintf('%s%s', scan_path, img_names_mat(round(num_proj_used/2), :));
         im_int = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo );
        [~, ht(1), dt(1)] = FilterPixelGPU( im_int, filt_pix_par );
        %[~, ht(1), dt(1)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo ), pixel_filter_threshold_proj, 0, pixel_filter_radius);
        filt_pix_par.threshold_hot  = median( ht );
        filt_pix_par.threshold_dark = median( dt );
    end
    
    % Preallocation
    roi_proj = zeros( numel( flat_corr_area1 ), numel( flat_corr_area2 ) , num_proj_used, 'single');
    
    % Lateral shift indices
    if isscalar( offset_shift )
        x0 = ones( 1, num_proj_used );
        x1 = im_shape_raw(1) * x0;
    else
        x0 = offset_shift';
        x1 = x0 + im_shape_raw(1) - max( x0 );        
    end
    
    % Preallocation
    im_shape_cropbin1 = floor( (x1(1) - x0(1) + 1) / raw_bin );
    proj = zeros( im_shape_cropbin1, im_shape_binned2, num_proj_used, 'single');
    fprintf( ' Allocated memory: %.2f GiB,', Bytes( proj, 3 ) )       
    projs_to_use = zeros( 1, size( proj,3), 'logical' );
    num_zeros = zeros( 1, num_proj_used );
    
    % Parallel loop over projections
    startS = ticBytes( gcp );
    parfor nn = 1:num_proj_used
        
        % Read projection
        filename = sprintf('%s%s', scan_path, img_names_mat(nn,:));
        im_int = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype, im_trafo );
        
        % Filter pixel
        im_int = FilterPixelGPU( im_int, filt_pix_par);
        
        % Dark field correction
        im_int = im_int - dark;
        
        % Correlation ROI
        roi_proj(:,:,nn) = im_int(flat_corr_area1,flat_corr_area2);

        % Remove lateral shift & Binning
        xx = x0(nn):x1(nn);
        im_float_binned = Binning( im_int(xx,:), raw_bin) / raw_bin^2;
        
        % Count zeros
        num_zeros(nn) = sum( im_float_binned(:) < 1 );
        
        % Reject image if any pixel is zero
        projs_to_use(nn) = ~boolean( num_zeros(nn)  );
        
        % Assign image to stack
        proj(:, :, nn) = im_float_binned;
    end
    toc_bytes.read_proj = tocBytes( gcp, startS );
    
    % Delete empty projections
    zz = ~projs_to_use;
    if sum( zz(:) )
        fprintf( ' Deleting empty projections' )
        proj(:,:,zz) = [];
    end
    if offset_shift ~= 0
        offset_shift(~projs_to_use) = [];
        x0(~projs_to_use) = [];
    end
    par.offset_shift_x0 = x0;
    par.offset_shift_x1 = x1;
    if vert_shift ~= 0
        vert_shift(~projs_to_use) = [];
    end
    tomo.vert_shift = vert_shift;
    if ~isempty( scan_position )
        scan_position(~projs_to_use) = [];
    end
    tomo.scan_position = scan_position;
    
    raw_min = min( proj(:) );
    raw_max = max( proj(:) );

    % Ring current normalization
    if ring_current_normalization(1)
        proj_ind_from_filenames = proj_nums;
        proj_ind_from_log = [cur.proj(proj_range).ind];
        if isequal( proj_ind_from_filenames,  proj_ind_from_log )
            proj_rc = [cur.proj(proj_range).val];
            proj_rcm = mean( proj_rc(:) );
            scale_factor = 100 ./ shiftdim( proj_rc(projs_to_use), -1 );
            proj = bsxfun( @times, proj, scale_factor );
            
            % Plot ring current
            if par.visual_output
                name = 'PETRA III beam current: Interpolation at image time stamps';
                if exist( 'hrc', 'var' ) && isvalid( hrc )
                    figure(hrc)
                else
                    hrc = figure( 'Name', name, 'WindowState', 'maximized');
                end
                subplot(1,1,1);
                plot( ref_nums, ref_rc(:), '.',proj_nums, proj_rc(:), '.' )
                axis tight
                xlabel( 'image no.' )
                ylabel( 'current / mA' )
                title( name )
                legend( sprintf( 'flats, mean: %.2f mA', ref_rcm), sprintf( 'projs, mean: %.2f mA', proj_rcm) )
                drawnow
                CheckAndMakePath( fig_path )
                saveas( hrc, sprintf( '%s%s.png', fig_path, regexprep( hrc.Name, '\ |:', '_') ) );
            end
        else
            fprintf('\n WARNING: projections not normalized by ring current. Names read from directory and log-file are inconsistent.\n')
        end
    end
    
    raw_min2 = min( proj(:) );
    raw_max2 = max( proj(:) );
    num_empty = sum( ~projs_to_use(:) );
    num_proj_used = num_proj_used - num_empty;
    fprintf( ' done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    fprintf( '\n crop left  min/max : %u %u', min( x0 ), max( x0 ) )
    fprintf( '\n crop right min/max : %u %u', min( x1 ), max( x1 ) )
    fprintf( '\n image shape cropbin1 : %u', im_shape_cropbin1 )
    PrintVerbose( num_empty, '\n discarded empty projections : %u, %.2f%%', num_empty, 100*num_empty/size(proj,3) )
    if sum( num_zeros )
        fprintf( '\n projections with zeros :' )
        % print #zeros if not all pixels are zero
        for nn = 1:numel(num_zeros)
            if num_zeros(nn) ~= 0
                if isequal( num_zeros(nn), numel_im_roi_binned )
                    fprintf( ' %u', nn )
                else
                    fprintf( ' %u:%u', nn, num_zeros(nn) )
                end
            end
        end
    end
    fprintf( '\n hot- / dark-pixel filter threshold : %f, %f', filt_pix_par.threshold_hot, filt_pix_par.threshold_dark )
    fprintf( '\n global min/max of projs after filtering and binning:  %6g %6g', raw_min, raw_max)
    fprintf( '\n global min/max of projs after dark-field correction and ring current normalization:  %6g %6g', raw_min2, raw_max2)
    
    %% Projection/flat field correlation and flat field correction %%%%%%%%
    par.raw_roi = raw_roi;
    par.proj_range = proj_range;
    par.ref_range = ref_range;
    [proj, corr, toc_bytes] = proj_flat_correlation( proj, flat, image_correlation, par, write, roi_proj, roi_flat, toc_bytes );
    %%%% STOP HERE TO CHECK FLATFIELD CORRELATION MAPPING %%%%%%%%%%%%%%%%%
    %%%% use 'proj_flat_sequ' to show results of the correlation
            
    proj_min0 = min( proj(:) );
    proj_max0 = max( proj(:) );
    fprintf( '\n global min/max after flat-field corrected:  %6g %6g', proj_min0, proj_max0);
    
    %% Plot data transfer from/to workers
    if par.visual_output 
        bytes_sum = [0 0];
        f = figure( 'Name', 'Parallel pool data transfer during image correlation', 'WindowState', 'maximized');
        Y = [];
        str = {};
        fn = fieldnames( toc_bytes );
        for nn = 1:numel( fn )
            fnn = fn{nn};
            X = toc_bytes.(fnn );
            Y = cat(2, Y, X );
            bytes_sum = bytes_sum + sum( X );
            fnn = regexprep( fnn, '_', ' ' );
            str = cat(2, str, { sprintf( '%s: to', fnn ), sprintf( '%s: from', fnn ) } );
        end
        plot( Y / 1024^3, 'o-' )
        axis tight
        title( sprintf( 'Data transfer of workers in parpool. Total: %.1f GiB (to), %.1f GiB (from)', bytes_sum / 1024^3 ) )
        xlabel( 'worker no.' )
        ylabel( 'transferred data / GiB' )
        legend( str )
        drawnow
        pause( 0.1 )
        CheckAndMakePath( fig_path )
        saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
    end
    
    %% Filter strong/full absorption (combine with iterative reco methods)
    if strong_abs_thresh < 1
        t = toc;
        fprintf( '\n set flat-corrected values below %f to one, ', strong_abs_thresh)
        parfor nn = 1:size( proj, 3 )
            im = proj(:,:,nn);
            m = im < strong_abs_thresh;
            im(m) = 0;
            im(m) = mean2( im );
            proj(:,:,nn) = im;
        end
        fprintf( ' done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    end
    
    fprintf( '\n sinogram size = [%g, %g, %g]', size( proj ) )
    
    %%% Angles %%
    if exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
        % for KIT cam this includes missing angles
        angles = [cur.proj.angle] / 180 * pi;
        if strcmpi(cam, 'kit')
            % drop angles where projections are missing
            angles = angles(1 + proj_nums);
        else
            angles = angles(proj_range);
        end
    elseif exist( h5log, 'file')
        angles = s_rot.value( ~boolean( stimg_key.value(logpar.n_dark+1:end) ) ) * pi / 180;
        angles = angles(proj_range);
    else
        num_proj = logpar.num_proj;
        switch lower( cam )
            case 'ehd'
                angles = tomo.rot_angle.full_range * (0:num_proj - 1) / (num_proj - 1); % EHD: ok
            case 'kit'
                angles = tomo.rot_angle.full_range * (0:num_proj - 1) / num_proj; % KIT: ok if logpar.projections exist
        end
    end
    % drop angles where projections are empty
    angles(~projs_to_use) = [];
    
    %% Figure & Save: sino slice
    nn = round( size( proj, 2) / 2);
    if isscalar( vert_shift )
        sino = squeeze( proj(:,nn,:) );
    else
        y = round( nn + vert_shift );
        y( y < 1 ) = 1;
        y( y > size( proj, 2) ) = size( proj, 2);
        for mm = size( proj, 3 ):-1:1
            sino(:,mm) = proj(:,y(mm),mm);
        end
    end
    
    CheckAndMakePath( reco_path )
    filename = sprintf( '%ssino_middle.tif', reco_path );
    write32bitTIFfromSingle( filename, sino);
    
    if par.visual_output
        if exist( 'h1' , 'var' ) && isvalid( h1 )
            figure(h1)
        else
            h1 = figure( 'Name', 'data and flat-and-dark-field correction', 'WindowState', 'maximized');
        end
        
        subplot(2,3,4)
        imsc1( FilterOutlier( proj(:,:,1), 0.005 ) )
        xticks([])
        yticks([])
        title(sprintf('intensity: first proj'))
        colorbar
        axis equal tight
        
        subplot(2,3,5)
        imsc1( FilterOutlier( proj(:,:,round(size(proj,3)/2)), 0.005 ) )
        xticks([])
        yticks([])
        title(sprintf('intensity: middle proj'))
        colorbar
        axis equal tight
        
        subplot(2,3,6)
        imsc1( FilterOutlier( proj(:,:,end), 0.005 ) )
        xticks([])
        yticks([])
        title(sprintf('intensity: last proj'))
        colorbar
        axis equal tight
        
        CheckAndMakePath( fig_path )
        saveas( h1, sprintf( '%s%s.png', fig_path, regexprep( h1.Name, '\ |:', '_') ) );
        
        f = figure( 'Name', 'sinogram', 'WindowState', 'maximized');        
        imsc1( sino )
        title( sprintf('sinogram: proj(:,%u,:)', nn) )
        colorbar
        axis equal tight
        
        saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
        drawnow
    end
    
    %% Ring artifact filter
    if ring_filter.apply && ring_filter.apply_before_stitching
        proj = p05_filter_ring_artefacts( ring_filter, proj, angles, par);
    end
    proj_min = min( proj(:) );
    proj_max = max( proj(:) );
    
    %% Write corrected projections
    if write.flatcor
        t = toc;
        fprintf( '\nSave flat-corrected projections.')
        CheckAndMakePath( flatcor_path, write.deleteFiles, write.beamtimeID )
        parfor nn = 1:size( proj, 3 )
            filename = sprintf('%sproj_%06u.tif', flatcor_path, nn );
            write32bitTIFfromSingle(filename, rot90( proj(:, :, nn) ) );
        end
        fprintf( ' done in %.1f (%.2f min)', toc-t, (toc-t)/60)
    end
else
    % Start parallel CPU pool %%%
    [poolobj, par.poolsize] = OpenParpool( par.poolsize, par.use_cluster, [beamtime_path filesep 'scratch_cc']);

    t = toc;
    %% Read sinogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if read_sino(1)
        
        % Sino path
        if ~isempty( read_sino_folder )
            sino_path = [scan_path read_sino_folder filesep];
        end
        
        % Sinogram file names
        data_struct = dir( [sino_path filesep '*.tif'] );
        if isempty( data_struct )
            error('\n Sinograms not found.')
        else
            sino_names = {data_struct.name};
            im_shape_binned2 = numel( sino_names );
        end
        sino_names_mat = NameCellToMat( sino_names );

        % Parameters
        filename = sprintf('%s%s', sino_path, sino_names_mat( round( im_shape_binned2 / 2 ), :));
        sino = read_image( filename );
        sino = read_sino_trafo( sino );
        im_shape_cropbin1 = size( sino, 2 );
        num_proj_read = size( sino, 1 );
        num_proj_used = num_proj_read;        

        % Angles
        if  ~exist( 'angles', 'var' )
            if isempty( tomo.rot_angle.full_range )
                cprintf( 'Red', '\nEnter full angle of rotation (including one additional increment) or vector of angles, in radians: ' );
                tomo.rot_angle.full_range = input( '' );
            end
            if isscalar( tomo.rot_angle.full_range )
                angles = tomo.rot_angle.full_range * (0:num_proj_read - 1) / num_proj_read;
            else
                angles = tomo.rot_angle.full_range;
            end
            if length( angles ) ~= num_proj_read
                error( 'Number of angles (%u) entered not consistent with sinogram (%u) read.', numel( angles), num_proj_read )
            end
        end
        
        % Preallocation
        fprintf( '\n Read sinograms.')
        proj = zeros( im_shape_cropbin1, im_shape_binned2, num_proj_read, 'single');
        fprintf( ' Allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
        
        % Read sinogram
        parfor nn = 1:size( proj, 2 )
            filename = sprintf('%s%s', sino_path, sino_names_mat(nn, :));
            sino = read_image( filename );
            sino = read_sino_trafo( sino );
            proj(:, nn, :) = permute( shiftdim( sino, -1 ) , [3 1 2] );
        end
        
        % Plot data
        if par.visual_output
            f = figure( 'Name', 'projection and sinogram', 'WindowState', 'maximized');
            
            subplot(1,2,1)
            imsc1( proj(:,:,1))
            title(sprintf('intensity: proj(:,:,1)'))
            colorbar
            axis equal tight
            
            subplot(1,2,2)
            nn = round( size( proj, 2) / 2);
            sino = squeeze( proj(:,nn,:) );
            sino = FilterOutlier( sino );
            imsc1( sino )
            title(sprintf('sinogram: proj(:,%u,:)', nn))
            colorbar
            axis equal tight
            
            drawnow
            
            CheckAndMakePath( fig_path )
            saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
            
        end
    end
    
    %% Read flat corrected projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if read_flatcor(1)
        
        if isempty( read_flatcor_path )
            read_flatcor_path = flatcor_path;
        end
        % File names
        data_struct = dir( [read_flatcor_path filesep 'proj*.*'] );
        if isempty( data_struct )
            fprintf('\n No flat corrected projections found! Switch to standard pre-processing.')
            read_flatcor = 0;
        else
            proj_names = {data_struct.name};
            num_proj_read = numel(proj_names);
            if num_proj_used > num_proj_read
                fprintf('\n Less projections available (%g) than demanded (%g)! Switch to standard pre-processing.', num_proj_read, num_proj_used )
                read_flatcor = 0;
            end
        end
        fprintf( '\n Read flat corrected projections.')
        if num_proj_read ~= num_proj_used
            fprintf('\n WARNING: Number of flat corrected projections read (%g) differs from number of projections to be processed (%g)!\n', num_proj_read, num_proj_used)
        end
        
        % Preallocation
        proj = zeros( im_shape_cropbin1, im_shape_binned2, num_proj_read, 'single');
        fprintf( ' Allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
        proj_names_mat = NameCellToMat( proj_names );
        
        % Read flat corrected projections
        parfor nn = 1:num_proj_read
            filename = sprintf('%s%s', flatcor_path, proj_names_mat(nn, :));
            im = read_image( filename );
            im = read_flatcor_trafo( im );
            proj(:, :, nn) = im;
        end
    end
    fprintf( ' done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
end

%% Phase retrieval before interactive mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tint_phase = 0;
if phase_retrieval.apply
    phase_retrieval.energy = energy;
    phase_retrieval.sample_detector_distance = sample_detector_distance;
    phase_retrieval.eff_pixel_size_binned = eff_pixel_size_binned;
    if isempty( tomo.take_neg_log )
        tomo.take_neg_log = 0;
    end
    if phase_retrieval.apply_before
        % Retrieval
        [proj, write, tint_phase] = phase_retrieval_func( proj, phase_retrieval, tomo, write, interactive_mode, par );
        tomo.rot_axis.position = tomo.rot_axis.position / phase_bin;
        tomo.rot_axis.offset = tomo.rot_axis.offset / phase_bin;
        [tomo.vol_shape, tomo.vol_size] = volshape_volsize( proj, tomo.vol_shape, tomo.vol_size, tomo.rot_axis.offset, verbose );
    end
end

%% TOMOGRAPHY: interactive mode to find rotation axis offset and tilt %%%%%
[tomo, angles, tint] = interactive_mode_rot_axis( par, logpar, phase_retrieval, tomo, write, interactive_mode, proj, angles);

%% Stitch projections
if par.stitch_projections
    t = toc;
    fprintf( '\nStitch projections:')
    % last projection within [0,pi)
    [~, num_proj_sti] = min( abs(angles - pi));
    % number of stitched projections
    num_proj_sti = num_proj_sti - 1;
    % index range of projections to be stitched
    xl = 1:round(tomo.rot_axis.position);
    xr = 1:xl(end)-1;
    im_shape_sti1 = numel( xl ) + numel( xr );
    % Preallocation
    proj_sti = zeros( im_shape_sti1 , size( proj, 2 ), num_proj_sti, 'single');
    for nn = 1:num_proj_sti
        nn2 = mod(num_proj_sti + nn - 1, size( proj, 3) ) + 1;
        im = zeros( im_shape_sti1, size( proj, 2 ));
        switch lower( par.stitch_method )
            case 'step'
                im = cat(1, proj(xl,:,nn), flipud( proj(xr,:,nn2) ) );
            case {'linear', 'sine'}
                % overlap region
                overlap = round(2 * tomo.rot_axis.position) - im_shape_cropbin1 : im_shape_cropbin1;
                % overlap ramp
                x = (0:1/(numel(overlap)-1):1);
                % 1D weight
                w = ones(im_shape_cropbin1, 1);
                switch lower( par.stitch_method )
                    case 'linear'
                        w(overlap) = 1 - x;
                    case 'sine'
                        w(overlap) = 0.5 * cos(pi*x) + 0.5;
                end
                % weighted projections
                iml = bsxfun(@times, proj(:,:,nn), w);
                imr = flipud( bsxfun(@times, proj(:,:,nn2), w ) );
                % stitched projection
                im(1:im_shape_cropbin1,:) = iml;
                im(end - im_shape_cropbin1 + 1:end,:) = im(end - im_shape_cropbin1 + 1:end,:) + imr;
        end
        proj_sti(:,:,nn) = im;
    end
    pause(0.01)
    proj = proj_sti;
    clear proj_sti;
    angles = angles(1:num_proj_sti);
    fprintf( ' done in %.1f (%.2f min)', toc-t, (toc-t)/60)
    fprintf( '\n shape of stitched projections : %u %u %u', size( proj ) )
    fprintf( '\n memory allocated : %.2f GiB', Bytes( proj, 3 ) )
end

%% Ring artifact filter %%
if ring_filter.apply && ~ring_filter.apply_before_stitching
    proj = p05_filter_ring_artefacts( ring_filter, proj, angles, par );
    proj_min = min( proj(:) );
    proj_max = max( proj(:) );
end

%% Crop projections at rotation axis position %%
if par.crop_at_rot_axis
    t = toc;
    fprintf( '\nCropping projections:')
    % Crop projections to avoid oversampling for scans with excentric rotation axis
    % and reconstruct WITHOUT stitching
    % Crop relative to rot axis position
    r = tomo.rot_axis.position / im_shape_cropbin1;
    if r < 0.5
        proj( 1:floor(tomo.rot_axis.position)-1, :, :) = [];
    else
        proj( ceil(tomo.rot_axis.position) + 1:end, :, :) = [];
    end
    if isempty( tomo.vol_shape )
        tomo.vol_shape = [raw_im_shape_binned1, raw_im_shape_binned1, raw_im_shape_binned2];
    end
    fprintf( ' done in %.1f (%.2f min)', toc-t, (toc-t)/60)
end

%% Save sinogram %%
if write.sino
    t = toc;
    fprintf( '\nSave sinogram:')
    CheckAndMakePath(sino_path, write.deleteFiles, write.beamtimeID)
    parfor nn = 1:size( proj, 2 )
        filename = sprintf( '%ssino_%06u.tif', sino_path, nn);
        sino = squeeze( proj( :, nn, :) )';
        write32bitTIFfromSingle( filename, sino )
    end
    pause(0.01)
    fprintf( ' done in %.1f s (%.2f min)', toc-t, (toc-t)/60)
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if phase_retrieval.apply
    if ~phase_retrieval.apply_before
        % Retrieval
        [proj, write, tint_phase] = phase_retrieval_func( proj, phase_retrieval, tomo, write, interactive_mode, par );
        % Post phase retrieval binning
        tomo.rot_axis.position = tomo.rot_axis.position / phase_bin;
        tomo.rot_axis.offset = tomo.rot_axis.offset / phase_bin;
        [tomo.vol_shape, tomo.vol_size] = volshape_volsize( proj, tomo.vol_shape, tomo.vol_size, tomo.rot_axis.offset, verbose );
    end
end

%% Tomographic reco %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tomo.run
    fprintf( '\nTomographic reconstruction:')
    fprintf( '\n method : %s', tomo.algorithm )
    fprintf( '\n angle scaling : %g', tomo.angle_scaling )
    fprintf( '\n angles [first last]/pi : [%g %g]', angles( [1 end] ) / pi )
    fprintf( '\n volume shape : [%g, %g, %g]', tomo.vol_shape )
    vol_mem = prod( tomo.vol_shape ) * 4;
    fprintf( '\n volume memory : %.2f GiB', vol_mem / 1024^3 )
    
    % Change 'reco_mode' to 'slice' if low on memory
    [mem_free, mem_avail, mem_total] = free_memory;
    fprintf( '\n system memory: free, available, total : %.3g GiB, %.3g GiB, %.3g GiB', mem_free/1024^3, mem_avail/1024^3, mem_total/1024^3)
    if vol_mem > 0.8 * mem_avail
        tomo.reco_mode = 'slice';
        fprintf( '\nSwitch to slice-wise reco (tomo.reco_mode = ''%s'') due to limited memory ( avail : %.1f GiB, vol : %.1f GiB) .', tomo.reco_mode, mem_avail / 1024^3, vol_mem / 1024^3 )
    end
    
    if par.stitch_projections
        rot_axis_offset_reco = 0;
    elseif par.crop_at_rot_axis
        rot_axis_offset_reco = tomo.rot_axis.position - size( proj, 1) / 2;
    else
        rot_axis_offset_reco = tomo.rot_axis.offset;
    end
    
    if isempty( tomo.take_neg_log )
        tomo.take_neg_log = 1;
    end
    
    % Delete redundant projection and angle
    if isequal( angles(1), angles(end) )
        angles(end) = [];
        proj(:,:,end) = [];
    end
    tomo.angles = tomo.rot_angle.offset + angles;
    
    % Filter sinogram
    if strcmpi( tomo.algorithm, 'fbp' )
        %%% Provide butterworth filtering also for iterative reconstruction !!!!
        fprintf( '\n Filter sino:' )
        t2 = toc;
        filt = iradonDesignFilter(tomo.fbp_filter.type, (1 + tomo.fbp_filter.padding) * size( proj, 1), tomo.fbp_filter.freq_cutoff);
        if tomo.butterworth_filter.apply(1)
            [b, a] = butter(tomo.butterworth_filter.order, tomo.butterworth_filter.frequ_cutoff);
            bw = freqz(b, a, numel(filt) );
            filt = filt .* bw;
        end
        proj_shape1 = size( proj, 1);
        take_neg_log = tomo.take_neg_log;
        padding = tomo.fbp_filter.padding;
        padding_method = tomo.fbp_filter.padding_method;
        parfor nn =  1:size( proj, 2)
            im = proj(:,nn,:);
            im = padarray( NegLog(im, take_neg_log), padding * [proj_shape1 0 0], padding_method, 'post' );
            im = real( ifft( bsxfun(@times, fft( im, [], 1), filt), [], 1, 'symmetric') );
            proj(:,nn,:) = im(1:proj_shape1,:,:);
        end
        pause(0.01)
        fprintf( ' done in %.2f min.', (toc - t2) / 60)
    else
        proj = NegLog( proj, tomo.take_neg_log );
    end
    
    % half weight pixel at rot axis pos as it is backprojected twice
    if par.crop_at_rot_axis
        r = tomo.rot_axis.position / im_shape_cropbin1;
        if r < 0.5
            proj( 1, :, :) = 0.5 * proj( 1, :, :) ;
        else
            proj( end, :, :) = 0.5 * proj( end, :, :) ;
        end
    end
    
    % Backprojection
    switch lower( tomo.reco_mode )
        case '3d'
            vol = zeros( tomo.vol_shape, 'single' );
            fprintf( '\n volume memory allocated for ''3D'' mode: %.2f GiB', Bytes( vol, 3 ) )
            
            fprintf( '\n Backproject:')
            t2 = toc;
            %%% Move to appropriate position or replace globally !!!!!!!!!!
            tomo.tilt_camera = ~interactive_mode.lamino * tomo.rot_axis.tilt;
            tomo.tilt_lamino = interactive_mode.lamino * tomo.rot_axis.tilt;
            %tomo.angles = tomo.rot_angle.offset + angles;
            tomo.rot_axis.offset = rot_axis_offset_reco;
            vol = astra_parallel3D( tomo, permute( proj, [1 3 2]) );
            pause(0.01)
            fprintf( ' done in %.2f min.', (toc - t2) / 60)
            
            vol_min = min( vol(:) );
            vol_max = max( vol(:) );
            
            %% Show orthogonal vol cuts
            if par.visual_output
                
                f = figure( 'Name', 'Volume cut z', 'WindowState', 'maximized');
                nn = round( size( vol, 3 ) / 2);
                im = squeeze( vol(:,:,nn) );
                im =  FilterOutlier( im, 0.01);
                imsc( im )
                axis equal tight
                title( sprintf( 'vol z = %u', nn ) )
                colorbar
                saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
                
                f = figure( 'Name', 'Volume cut y', 'WindowState', 'maximized');
                nn = round( size( vol, 2 ) / 2);
                im = rot90( squeeze( vol(:,nn,:) ), -2);
                im = FilterOutlier( im, 0.01 );
                if size( im , 1) < size( im , 2)
                    imsc( im )
                    title( sprintf( 'vol y = %u', nn ) )
                else
                    imsc( im' )
                    title( sprintf( 'vol y = %u, rotated', nn ) )
                end
                axis equal tight
                colorbar
                saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
                
                f = figure( 'Name', 'Volume cut x', 'WindowState', 'maximized');
                nn = round( size( vol, 1 ) / 2);
                im = flipud( squeeze( vol(nn,:,:) ) );
                im =  FilterOutlier( im, 0.01);
                if size( im , 1) < size( im , 2)
                    imsc( im )
                    title( sprintf( 'vol x = %u', nn ) )
                else
                    imsc( im' )
                    title( sprintf( 'vol x = %u rotated', nn ) )
                end
                axis equal tight
                colorbar
                saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
                
                drawnow
            end
            
            if phase_retrieval.apply
                reco_path = write.reco_phase_path;
            end
            CheckAndMakePath( reco_path, 0 )
            
            %% Save ortho slices
            
            % Save ortho slices x
            nn = round( size( vol, 1 ) / 2);
            im = squeeze( vol(nn,:,:) );
            filename = sprintf( '%sreco_xMid.tif', reco_path );
            write32bitTIFfromSingle( filename, rot90(im,-1) );
            
            % Save ortho slices y
            nn = round( size( vol, 2 ) / 2);
            im = squeeze( vol(:,nn,:) );
            filename = sprintf( '%sreco_yMid.tif', reco_path );
            write32bitTIFfromSingle( filename, rot90(im,-1) );
            
            % Save ortho slices z
            nn = round( size( vol, 3 ) / 2);
            im = squeeze( vol(:,:,nn) );
            filename = sprintf( '%sreco_zMid.tif', reco_path );
            write32bitTIFfromSingle( filename, rot90(im,0) );
            
            % Save volume
            if ~isfield( write, 'reco_binning_factor')
                write.reco_binning_factor = 2;
            end
            if write.reco
                reco_bin = write.reco_binning_factor; % alias for readablity
                CheckAndMakePath( reco_path, 0 )
                
                % Save reco path to file
                filename = [userpath, filesep, 'experiments/p05/path_to_latest_reco'];
                fid = fopen( filename , 'w' );
                fprintf( fid, '%s', reco_path );
                fclose( fid );
                
                % Single precision: 32-bit float tiff
                write_volume( write.float, vol, 'float', reco_path, raw_bin, phase_bin, 1, 0, verbose, '', write.deleteFiles, write.beamtimeID);
                
                % Compression of dynamic range
                if write.uint8 || write.uint8_binned || write.uint16 || write.uint16_binned
                    [tlow, thigh] = compression( vol, write.compression.method, write.compression.parameter );
                else
                    tlow = 0;
                    thigh = 1;
                end
                
                % 16-bit tiff
                write_volume( write.uint16, (vol - tlow)/(thigh - tlow), 'uint16', reco_path, raw_bin, phase_bin, 1, 0, verbose, '', write.deleteFiles, write.beamtimeID);
                
                % 8-bit tiff
                write_volume( write.uint8, (vol - tlow)/(thigh - tlow), 'uint8', reco_path, raw_bin, phase_bin, 1, 0, verbose, '', write.deleteFiles, write.beamtimeID);
                
                % Bin data
                if write.float_binned || write.uint16_binned || write.uint8_binned || write.uint8_segmented
                    fprintf( '\n Binning:')
                    t2 = toc;
                    vol = Binning( vol, reco_bin ) / reco_bin^3;
                    fprintf( ' done in %.2f min.', (toc - t2) / 60)
                end
                
                % Binned single precision: 32-bit float tiff
                write_volume( write.float_binned, vol, 'float', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose, '', write.deleteFiles, write.beamtimeID);
                
                % 16-bit tiff binned
                write_volume( write.uint16_binned, (vol - tlow)/(thigh - tlow), 'uint16', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose, '', write.deleteFiles, write.beamtimeID);
                
                % 8-bit tiff binned
                write_volume( write.uint8_binned, (vol - tlow)/(thigh - tlow), 'uint8', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose, '', write.deleteFiles, write.beamtimeID);
                
                % segmentation
                if write.uint8_segmented
                    [vol, out] = segment_volume(vol, 2^10, par.visual_output, verbose);
                    save_path = write_volume( 1, vol/255, 'uint8', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose, '_segmented', write.deleteFiles, write.beamtimeID);
                    save( sprintf( '%ssegmentation_info.m', save_path), 'out', '-mat', '-v7.3')
                end
            end
            
        case {'slice', '2d'}
            %% Slicewise backprojection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
              % Remove vertical shift for spiral CT
            if vert_shift
                t2 = toc;
                fprintf( '\n Remove vertical shift:' )
                parfor nn = 1:size( proj, 3 )
                    im = proj(:,:,nn);
                    imt = imtranslate( im, [-vert_shift(nn) 0], 'linear' );
                    proj(:,:,nn) = imt;
                end
                fprintf( ' done in %.1f s (%.2f min)', toc-t2, (toc-t2)/60 )
            end
            
            fprintf( '\n Backproject and save slices:')
            t2 = toc;
            
            if phase_retrieval.apply
                reco_path = write.reco_phase_path;
            end
            
            if write.reco
                reco_bin = write.reco_binning_factor; % alias for readablity
                CheckAndMakePath( reco_path, 0 )
                % Save reco path to file
                filename = [userpath, filesep, 'experiments/p05/path_to_latest_reco'];
                fid = fopen( filename , 'w' );
                fprintf( fid, '%s', reco_path );
                fclose( fid );
            end
            
            % tomo.angles = tomo.rot_angle.offset + angles;
            tomo.rot_axis.offset = rot_axis_offset_reco;
            
            % Reconstruct central slices first
            [~, ind] = sort( abs( (1:size( proj, 2)) - round( size( proj, 2) / 2 ) ) );
            
            % Loop over slices
            vol_min = Inf;
            vol_max = -Inf;
            indshow = [1, round( [0.25 0.5 0.75 1] * size( proj, 2) )];
            for nn = 1:size( proj, 2 )
                
                % Backproject
                vol = rot90( astra_parallel2D( tomo, permute( proj(:,ind(nn),:), [3 1 2]) ), -1);
                
                vol_min = min( vol_min, min( vol(:) ) );
                vol_max = max( vol_max, max( vol(:) ) );
                
                % Show orthogonal vol cuts
                if par.visual_output
                    if sum( ind(nn) == indshow )
                        figure( 'Name', 'Volume slice', 'WindowState', 'maximized');
                        imsc( vol )
                        axis equal tight
                        title( sprintf( 'vol z = %u', ind(nn) ) )
                        colorbar
                        drawnow
                        pause( 0.1 )
                    end
                end
                
                % Save ortho slices z
                if nn == round( tomo.vol_shape / 2 )
                    filename = sprintf( '%sreco_zMid.tif', reco_path );
                    write32bitTIFfromSingle( filename, rot90(vol,0) );
                end
                
                % Save
                %%% ADD other format options!!!!!!!!!!!!!!!!!!!!!!!!!
                if write.reco
                    % Single precision: 32-bit float tiff
                    write_volume( write.float, vol, 'float', reco_path, raw_bin, phase_bin, 1, ind(nn) - 1, 0, '', write.deleteFiles, write.beamtimeID);
                end
            end
    end
    fprintf( ' \n tomo reco done in %.1f s (%.2f min)', toc-t, (toc-t)/60 )
end

%% Write reco log file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CheckAndMakePath( reco_path )
save( sprintf( '%sangles.mat', reco_path), 'angles' );
save( sprintf( '%soffset_shift.mat', reco_path), 'offset_shift' );
if write.reco
    logfile_path = reco_path;
    
else
    logfile_path = write.path;
end
CheckAndMakePath( logfile_path )
if write.reco
    if phase_retrieval.apply
        logfile_name = sprintf( '%sreco_phase_%s_rawBin%u.log', logfile_path, write.phase_appendix, raw_bin );
    else
        logfile_name = sprintf( '%sreco_rawBin%u.log', logfile_path, raw_bin);
    end
    fid = fopen( logfile_name, 'w');
    fprintf(fid, 'scan_name : %s\n', scan_name);
    fprintf(fid, 'beamtime_id : %s\n', beamtime_id);
    fprintf(fid, 'scan_path : %s\n', scan_path);
    fprintf(fid, 'reco_path : %s\n', reco_path);
    fprintf(fid, 'MATLAB notation, index of first element: 1, range: first:stride:last\n');
    fprintf(fid, 'MATLAB version : %s\n', version);
    fprintf(fid, 'Git commit ID : %s\n', git_commit_id);
    fprintf(fid, 'platform : %s\n', computer);
    
    fprintf(fid, 'effective_pixel_size : %g micron\n', eff_pixel_size * 1e6);
    fprintf(fid, 'effective_pixel_size_binned : %g micron\n', eff_pixel_size_binned * 1e6);
    fprintf(fid, 'energy : %g eV\n', energy);
    fprintf(fid, 'sample_detector_distance : %f m\n', sample_detector_distance);
    if ~read_sino
        fprintf(fid, 'camera : %s\n', cam);
        fprintf(fid, 'exposure_time : %f\n', exposure_time);
        fprintf(fid, 'num_dark_found : %u\n', num_dark);
        fprintf(fid, 'num_ref_found : %u\n', num_ref_found);
        fprintf(fid, 'num_ref_used : %u\n', num_ref_used);
        fprintf(fid, 'ref_range : %u:%u:%u\n', ref_range(1), ref_range(2) - ref_range(1), ref_range(end) );
        fprintf(fid, 'num_proj_found : %u\n', num_proj_found);
        fprintf(fid, 'num_proj_used : %u\n', num_proj_used);
        fprintf(fid, 'proj_range : %u:%u:%u\n', proj_range(1), proj_range(2) - proj_range(1), proj_range(end) );
        fprintf(fid, 'im_shape_raw : %u %u\n', im_shape_raw);
        fprintf(fid, 'raw_roi : ');
        fprintf(fid, ' %u', raw_roi);
        fprintf(fid, '\n' );
        fprintf(fid, 'im_shape_roi : %u %u\n', size( im_roi ));
        fprintf(fid, 'im_shape_binned : %u %u\n', im_shape_binned1, im_shape_binned2);
        fprintf(fid, 'im_shape_cropbin1 : %u %u\n', im_shape_cropbin1);
        fprintf(fid, 'raw_binning_factor : %u\n', raw_bin);
        fprintf(fid, 'ring_current_normalization : %u\n', ring_current_normalization);
        fprintf(fid, 'image_correlation.method : %s\n', image_correlation.method);
        fprintf(fid, 'image_correlation.num_flats : %u\n', image_correlation.num_flats);
        fprintf(fid, 'image_correlation.area_width : %u:%u:%u\n', image_correlation.area_width(1), image_correlation.area_width(2) - image_correlation.area_width(1), image_correlation.area_width(end));
        fprintf(fid, 'image_correlation.area_height : %u:%u:%u\n', image_correlation.area_height(1), image_correlation.area_height(2) - image_correlation.area_height(1), image_correlation.area_height(end));
    end
    if ~read_flatcor && ~read_sino
        fprintf(fid, 'min_max_of_all_darks : %6g %6g\n', darks_min, darks_max);
        fprintf(fid, 'min_max_of_median_dark : %6g %6g\n', dark_med_min, dark_med_max);
        fprintf(fid, 'min_max_of_all_flats : %6g %6g\n', flat_min, flat_max);
        fprintf(fid, 'min_max_of_all_corrected_flats : %6g %6g\n', flat_min2, flat_max2);
        fprintf(fid, 'min_max_of_all_raws :  %6g %6g\n', raw_min, raw_max);
        fprintf(fid, 'min_max_of_all_corrected_raws :  %6g %6g\n', raw_min2, raw_max2);
        fprintf(fid, 'min_max_of_all_flat_corr_projs : %g %g \n', proj_min, proj_max);
    end
    fprintf(fid, 'strong_abs_thresh : %f m\n', strong_abs_thresh);
    % Phase retrieval
    fprintf(fid, 'phase_retrieval.apply : %u\n', phase_retrieval.apply);
    if phase_retrieval.apply
        fprintf(fid, 'phase_retrieval.padding : %u\n', phase_retrieval.padding);
        fprintf(fid, 'phase_retrieval.method : %s\n', phase_retrieval.method);
        fprintf(fid, 'phase_retrieval.regularisation_parameter : %f\n', phase_retrieval.reg_par);
        fprintf(fid, 'phase_retrieval.binary_filter_threshold : %f\n', phase_retrieval.bin_filt);
        fprintf(fid, 'phase_retrieval.cutoff_frequency : %f pi\n', phase_retrieval.cutoff_frequ / pi);
        fprintf(fid, 'phase_retrieval.post_binning_factor : %u\n', phase_bin);
    end
    if tomo.run
        % Volume
        fprintf(fid, 'tomo.vol_shape : %u %u %u\n', tomo.vol_shape(1), tomo.vol_shape(2), tomo.vol_shape(3));
        fprintf(fid, 'tomo.vol_size : %f %f %f %f %f %f\n', tomo.vol_size(1), tomo.vol_size(2), tomo.vol_size(3), tomo.vol_size(4), tomo.vol_size(5), tomo.vol_size(6));
        % Rotation
        fprintf(fid, 'crop_at_rot_axis : %u\n', par.crop_at_rot_axis);
        fprintf(fid, 'par.stitch_projections : %u\n', par.stitch_projections);
        fprintf(fid, 'par.stitch_method : %s\n', par.stitch_method );
        fprintf(fid, 'tomo.reco_mode : %s\n', tomo.reco_mode);
        fprintf(fid, 'tomo.rot_angle.full_range : %f * pi rad\n', tomo.rot_angle.full_range / pi);
        fprintf(fid, 'tomo.rot_angle.offset : %f * pi rad\n', tomo.rot_angle.offset / pi);
        fprintf(fid, 'rot_axis_offset_reco : %f\n', rot_axis_offset_reco);
        fprintf(fid, 'tomo.rot_axis.position : %f\n', tomo.rot_axis.position);
        fprintf(fid, 'tomo.rot_axis.tilt : %f\n', tomo.rot_axis.tilt);
        fprintf(fid, 'raw_image_binned_center : %f\n', im_shape_cropbin1 / 2);
        fprintf(fid, 'interactive_mode.rot_axis_pos : %u\n', interactive_mode.rot_axis_pos);
        fprintf(fid, 'interactive_mode.phase_retrieval : %u\n', interactive_mode.phase_retrieval);
        % Ring filter
        fprintf(fid, 'ring_filter.apply : %u\n', ring_filter.apply);
        if ring_filter.apply
            fprintf(fid, 'ring_filter.method : %s\n', ring_filter.method);
            fprintf(fid, 'ring_filter.apply_before_stitching : %u\n', ring_filter.apply_before_stitching);
            switch lower( ring_filter.method )
                case 'wavelet-fft'
                    fprintf(fid, 'ring_filter.waveletfft.wname : %s\n', ring_filter.waveletfft.wname );
                    fprintf(fid, 'ring_filter.waveletfft.dec_levels : [%s]\n', sprintf( ' %u', ring_filter.waveletfft.dec_levels) );
                    fprintf(fid, 'ring_filter.waveletfft.sigma : %f\n', ring_filter.waveletfft.sigma );
                case 'jm'
                    fprintf(fid, 'ring_filter.jm.median_width : %s\n', sprintf( '%u ', ring_filter.jm.median_width) );
            end
        end
        % Tomo
        fprintf(fid, 'tomo.algorithm : %s\n', tomo.algorithm );
        switch tomo.algorithm
            case 'fbp' %'fbp-astra'}
                fprintf(fid, 'tomo.fbp_filter.type : %s\n', tomo.fbp_filter.type);
                fprintf(fid, 'tomo.fbp_filter.freq_cutoff : %f\n', tomo.fbp_filter.freq_cutoff);
                fprintf(fid, 'tomo.fbp_filter.padding : %u\n', tomo.fbp_filter.padding);
                fprintf(fid, 'tomo.fbp_filter.padding_method : %s\n', tomo.fbp_filter.padding_method);
            case 'sirt'
                fprintf(fid, 'tomo.iterations : %u\n',  tomo.iterations );
                fprintf(fid, 'tomo.sirt.MinConstraint : %f\n',  tomo.sirt.MinConstraint);
                fprintf(fid, 'tomo.sirt.MaxConstraint : %f\n',  tomo.sirt.MaxConstraint);
            case 'cgls';'sart';'em';
                fprintf(fid, 'tomo.iterations : %u\n',  tomo.iterations );
        end
        fprintf(fid, 'tomo.butterworth_filter.apply : %u\n', tomo.butterworth_filter.apply);
        fprintf(fid, 'tomo.butterworth_filter.order : %u\n', tomo.butterworth_filter.order);
        fprintf(fid, 'tomo.butterworth_filter.frequ_cutoff : %f\n', tomo.butterworth_filter.frequ_cutoff);
        fprintf(fid, 'tomo.astra_pixel_size : %f\n', tomo.astra_pixel_size);
        fprintf(fid, 'tomo.take_neg_log : %u\n', tomo.take_neg_log);
        fprintf(fid, 'gpu_name : %s\n', gpu.Name);
        if exist( 'vol_min', 'var')
            fprintf(fid, '[volume_min volume_max] : [%g %g]\n', vol_min, vol_max);
        end
        fprintf(fid, 'write.float : %u\n', write.float);
        if strcmpi( tomo.reco_mode, '3d' )
            fprintf(fid, 'write.float_binned : %u\n', write.float_binned);
            fprintf(fid, 'write.uint16 : %g\n', write.uint16);
            fprintf(fid, 'write.uint16_binned : %g\n', write.uint16_binned);
            fprintf(fid, 'write.uint8 : %u\n', write.uint8);
            fprintf(fid, 'write.uint8_binned : %u\n', write.uint8_binned);
            if write.uint16 || write.uint8 || write.uint16_binned || write.uint8_binned
                fprintf(fid, 'write.compression.method : %s\n', write.compression.method);
            end
            if exist( 'tlow', 'var' ) && exist( 'thigh', 'var' )
                fprintf(fid, 'compression.limits : %f %f\n', tlow, thigh);
            end
            fprintf(fid, 'reco_bin : %u\n', write.reco_binning_factor);
        end
        fprintf(fid, 'full_reconstruction_time : %.1f s\n', toc);
        fprintf(fid, 'date_of_reconstruction : %s\n', datetime);
        fprintf(fid, 'tomo.rot_axis.offset at %u x binning : %f\n', raw_bin, rot_axis_offset_reco);
    end
    fclose(fid);
    % End of log file
    fprintf( '\n log file : \n%s', logfile_name)
    fprintf( '\n reco_path : \n%s', reco_path)
end
PrintVerbose( interactive_mode.rot_axis_pos, '\nTime elapsed in interactive rotation axis centering mode: %g s (%.2f min)', tint, tint / 60 );
PrintVerbose( interactive_mode.phase_retrieval, '\nTime elapsed in interactive phase retrieval mode: %g s (%.2f min)', tint_phase, tint_phase / 60 );
fprintf( '\nTime elapsed for computation: %g s (%.2f min)', toc - tint -tint_phase, (toc - tint - tint_phase) / 60 );
fprintf( '\nFINISHED: %s at %s\n', scan_name, datetime )
if exist('fast_reco','var') && fast_reco.run
    cprintf( 'Red', '\nATTENTION: fast reco mode was turned on!\n' )
end
fprintf( '\n')
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dbclear if error
