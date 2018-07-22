% P05 reconstruction pipeline: preprocessing, filtering, phase retrieval,
% tomographic reconstruction, simple segmentation, etc.
%
% USAGE:
% Set parameters in PARAMETERS / SETTINGS section below and run script.
%
% HOW TO RUN THE SCRIPT:
% - Editor windows: press 'F5' when focus is in the Editor window
% - Editor tab: click 'Run' in the
% - Command Window: type 'p05_reco' and Enter
%
% HOW TO AUTOMATICALLY LOOP OVER RECONSTRUCTIONS:
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
% Written by Julian Moosmann. First version: 2016-09-28. Last modifcation:
% 2018-06-21

if ~exist( 'external_parameter' ,'var')
    clearvars
end
close all hidden % close all open windows
%dbstop if error

%% PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fast_reco = 1; % OVERWRITES PARAMETERS SET BELOW
stop_after_data_reading(1) = 0; % for data analysis, before flat field correlation

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scan_path = ...
    '/asap3/petra3/gpfs/p05/2017/data/11004016/raw/syn002_6L_PEEK_4w_002';
read_flatcor = 0; % read flatfield-corrected images from disc, skips preprocessing
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = [1200 1500]; % if []: use full image; if [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2). When raw_roi(2) < 0 reads until end - |raw_roi(2)|; if negative scalar: auto roi, selects ROI automatically.Not working for *.raw data where images are flipped.
raw_bin = 4; % projection binning factor: integer
bin_before_filtering(1) = 1; % Apply binning before filtering pixel. less effective, but much faster especially for KIT camera.
excentric_rot_axis = 0; % off-centered rotation axis increasing FOV. -1: left, 0: centeerd, 1: right. influences tomo.rot_axis.corr_area1
crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
stitch_projections = 0; % for 2 pi scans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixe
ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';
% Correlation of projections and flat fields important for DCM data. For
% DMM it does not work well. Available methods:
% 'ssim-ml' : Matlab's structural similarity index (SSIM), includes Gaussian smoothing
% 'ssim' : own implementation of SSIM, smoothing not yet implemented, usually worse than 'ssim-ml'
% 'entropy' : entropy measure of proj over flat
% 'cov' : cross covariance
% 'corr' : cross correlation = normalized cross covariance
% 'std' : standard deviation of proj over flat
% 'diff': difference of proj and flat
% 'shift': computes relative shift from peak of cross-correlation map
% 'none' : no correlation, use median flat
% 'cross-entropy*' : variants of (asymmetric) cross entropy
image_correlation.num_flats = 1; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
image_correlation.area_width = [0 0.02];%[0.98 1];% correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.area_height = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.shift.max_pixelshift = 0.25; % maximum pixelshift allowed for 'shift'-correlation method: if 0 use the best match (i.e. the one with the least shift), if > 0 uses all flats with shifts smaller than image_correlation.shift.max_pixelshift
ring_filter.apply = 0; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft.dec_levels = 2:5; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft.wname = 'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter.jm.median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
strong_abs_thresh = 1; % Experimental: if 1: does nothing, if < 1: flat-corrected values below threshold are set to one
decimal_round_precision = 2; % precision when rounding pixel shifts
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 1.5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.15; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0;
tomo.reco_mode = '3D'; 'slice'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = []; %[-0.5 0.5 -0.5 0.5 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle.full_range = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle.offset = pi; % global rotation of reconstructed volume
tomo.rot_axis.offset = [];%-2.5;[];% if empty use automatic computation
tomo.rot_axis.position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis.tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
tomo.rot_axis.corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
tomo.rot_axis.corr_area2 = []; % ROI to correlate projections at angles 0 & pi
tomo.rot_axis.corr_gradient = 0; % use gradient of intensity maps if signal variations are too weak to correlate projections
tomo.fbp_filter.type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
tomo.fbp_filter.freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
tomo.fbp_filter.padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
tomo.fbp_filter.padding_method = 'symmetric';
tomo.butterworth_filter.apply = 0; % use butterworth filter in addition to FBP filter
tomo.butterworth_filter.order = 1;
tomo.butterworth_filter.frequ_cutoff = 0.9;
tomo.astra_pixel_size = 1; % detector pixel size for reconstruction: if different from one 'tomo.vol_size' must to be ajusted, too!
tomo.take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
tomo.algorithm = 'fbp';'sirt'; 'cgls';
tomo.iterations = 40; % for 'sirt' or 'cgls'.
tomo.sirt.MinConstraint = []; % If specified, all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions, for example.
tomo.sirt.MaxConstraint = []; % If specified, all values above MaxConstraint will be set to MaxConstraint.
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = '';'/gpfs/petra3/scratch/moosmanj';% absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 1; % write to 'scratch_cc' instead of 'processed'
write.parfolder = '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder.flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder.phase_map = ''; % subfolder in 'phase_map'
write.subfolder.sino = ''; % subfolder in 'sino'
write.subfolder.reco = ''; % subfolder in 'reco'
write.flatcor = 0; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % additionally save 16bit unsigned integer tiff using 'write.compression.method'
write.uint8 = 0; % additionally save binned 8bit unsigned integer tiff using 'write.compression.method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 0; % additionally save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % additionally save binned 16bit unsigned integer tiff using 'write.compression.method'
write.uint8_binned = 0; % additionally save binned 8bit unsigned integer tiff using 'wwrite.compression.method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression.method = 'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression.parameter = [0.20 0.15]; % compression-method specific parameter
% dynamic range is compressed s.t. new dynamic range assumes
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression.parameter, eg. [-0.01 1]
% 'std' : NUM = write.compression.parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression.parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 1; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
poolsize = 0.60; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
%%% EXPERIMENTAL OR NOT YET IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.uint8_segmented = 0; % experimental: threshold segmentation for histograms with 2 distinct peaks: __/\_/\__

%%% END OF PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reco loop: parameters set by external call from 'p05_reco_loop' %%%%%%%
if exist( 'external_parameter' ,'var')
    %verbose = 0;
    visual_output = 0;
    dbstop if error
    interactive_mode.rot_axis_pos = 0;
    fast_reco = 0;
    stop_after_data_reading(1) = 0;
    fn = fieldnames( external_parameter );
    for nn = 1:numel( fn )
        var_name = fn{nn};
        var_val = getfield( external_parameter, var_name );
        assignin('caller', var_name, var_val )
    end
    clear external_parameter;
end

% overwrite parameters for fast reconstruction
if exist('fast_reco','var') && fast_reco(1)
    raw_roi = [1000 -1000]; 
    raw_bin = 4; 
    bin_before_filtering(1) = 1;
    proj_range = 3;
    ref_range = 10;
    image_correlation.method = 'none';
    %tomo.reco_mode = 'slice';
    write.to_scratch = 1;
end

%% Preprocessing up to proj/flat correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if interactive_mode.rot_axis_tilt && strcmpi( tomo.reco_mode, 'slice' )
    error( 'Slicewise reconstruction and reconstruction with tilted rotation axis are not compatible!' )
end
% set (relative) vol_shape from vol_size if vol_shape is empty, but vol_size is given
if ~isempty( tomo.vol_size ) && isempty( tomo.vol_shape )
    tomo.vol_shape = tomo.vol_size(2:2:end) - tomo.vol_size(1:2:end);
end
imsc1 = @(im) imsc( rot90( im ) );
prnt = @(varargin) PrintVerbose( verbose, varargin{:});
prnt( 'Start reconstruction of ')
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';
astra_clear % if reco was aborted, ASTRA memory is not cleared
if ~isempty( tomo.rot_axis.offset ) && ~isempty( tomo.rot_axis.position )
    error('tomo.rot_axis.offset (%f) and tomo.rot_axis.position (%f) cannot be used simultaneously. One must be empty.', tomo.rot_axis.offset, tomo.rot_axis.position)
end
if abs(excentric_rot_axis)
    if crop_at_rot_axis(1) && stitch_projections(1)
        error( 'Do not use ''stitch projections'' in combination with ''crop_at_rot_axis''. Cropping at rot axis only makes sense without stitching in order to avoid oversampling artefacts within the overlap region.' )
    end
end

%%%% Folders: input %%%
while scan_path(end) == filesep
    scan_path(end) = [];
end
[raw_path, scan_name] = fileparts(scan_path);
% Save raw path to file for shell short cut
filename = [userpath, filesep, 'experiments/p05/path_to_latest_raw'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', raw_path );
fclose( fid );
% Scan path
scan_path = [scan_path, filesep];
[beamtime_path, raw_folder] = fileparts(raw_path);
[~, beamtime_id] = fileparts(beamtime_path);
if ~strcmp(raw_folder, 'raw')
    error('Given path does not contain a ''raw'' folder: %s', raw_folder)
end
prnt( '%s', scan_name)
prnt( '\n scan_path:%s', scan_path)
% Save scan path to file
filename = [userpath, filesep, 'experiments/p05/path_to_latest_scan'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', scan_path );
fclose( fid );

%%% Folders: output %%%
out_folder = 'processed';
if write.to_scratch
    out_folder = 'scratch_cc';
end
if isempty( write.path )
    write.path = [beamtime_path, filesep, out_folder, filesep, scan_name];
else
    write.path = [write.path, filesep, scan_name];
end
if ~isempty(write.parfolder)
    write.path = [write.path, filesep, write.parfolder];
end

% Path to flat-field corrected projections
if isempty( write.subfolder.flatcor )
    flatcor_path = [write.path, filesep, 'flat_corrected', filesep];
else
    flatcor_path = [write.path, filesep, 'flat_corrected', filesep, write.subfolder.flatcor, filesep];
end
PrintVerbose(verbose & write.flatcor, '\n flatcor_path: %s', flatcor_path)

% Path to retreived phase maps
if isempty( write.subfolder.phase_map )
    phase_map_path = [write.path, filesep, 'phase_map', filesep];
else
    phase_map_path = [write.path, filesep, 'phase_map', filesep, write.subfolder.phase_map, filesep];
end
PrintVerbose(verbose & write.phase_map, '\n phase_map_path: %s', phase_map_path)

% Sinogram path
if isempty( write.subfolder.sino )
    sino_path = [write.path, filesep, 'sino', filesep];
    sino_phase_path = [write.path, filesep, 'sino_phase', filesep];
else
    sino_path = [write.path, filesep, 'sino', filesep, write.subfolder.sino, filesep];
    sino_phase_path = [write.path, filesep, 'sino_phase', filesep, write.subfolder.sino, filesep];
end
PrintVerbose(verbose & write.sino, '\n sino_path: %s', sino_path)
PrintVerbose(verbose & write.phase_sino, '\n sino_phase_path: %s', sino_phase_path)

% Reco path
if isempty( write.subfolder.reco )
    reco_path = [write.path, filesep, 'reco', filesep];
else
    reco_path = [write.path, filesep, 'reco', filesep, write.subfolder.reco, filesep];
end

% Memory
prnt( '\n hostname : %s', getenv( 'HOSTNAME' ) );
[mem_free, mem_avail, mem_total] = free_memory;
prnt( '\n system memory: free, available, total : %.3g GiB, %.3g GiB, %.3g GiB', mem_free/1024^3, mem_avail/1024^3, mem_total/1024^3)
if isempty( tomo.astra_gpu_index )
    tomo.astra_gpu_index = 1:gpuDeviceCount;
end
for nn = tomo.astra_gpu_index
    gpu = parallel.gpu.GPUDevice.getDevice( nn );
    prnt( '\n gpu %u memory : available, total, percent = %.3g GiB, %.3g GiB, %.2f%%', nn, gpu.AvailableMemory/1024^3, gpu.TotalMemory/1024^3, 100*gpu.AvailableMemory/gpu.TotalMemory)
end

%%% File names %%%

% Projection file names
proj_names = FilenameCell( [scan_path, '*.img'] );
raw_data = 0;
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*proj*.tif'] );
    raw_data = 0;
end
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*img*.tif'] );
    raw_data = 0;
end
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*img*.raw'] );
    raw_data = 1;
end
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*proj*.raw'] );
    raw_data = 1;
end
num_proj_found = numel(proj_names);

% Ref file names
ref_names = FilenameCell( [scan_path, '*.ref'] );
if isempty( ref_names )
    ref_names = FilenameCell( [scan_path, '*ref*.tif'] );
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
num_ref_found = numel(ref_names);
if isempty( ref_range )
    ref_range = 1;
end
if numel( ref_range ) == 1
    ref_range = 1:ref_range:num_ref_found;
end
% position of running index
if strcmpi( ref_names{1}(end-6:end-4), 'ref' )
    imtype_str_flag = 0;
elseif strcmpi( ref_names{1}(end-11:end-9), 'ref' )
    imtype_str_flag = 1;
else
    imtype_str_flag = -1;
end
num_ref_used = numel( ref_range );
ref_names_mat = NameCellToMat( ref_names(ref_range) );
ref_nums = CellString2Vec( ref_names(ref_range) , imtype_str_flag);
prnt( '\n number of refs found : %g', num_ref_found)
prnt( '\n number of refs used : %g', num_ref_used)
prnt( '\n reference range used : %g:%g:%g%', ref_range(1), ref_range(2) - ref_range(1), ref_range(end))

% Dark file names
dark_names = FilenameCell( [scan_path, '*.dar'] );
if isempty( dark_names )
    dark_names = FilenameCell( [scan_path, '*dar*.tif'] );
end
if isempty( dark_names )
    dark_names =  FilenameCell( [scan_path, '*dar*.raw'] );
end
dark_nums = CellString2Vec( dark_names, imtype_str_flag );
num_dark = numel(dark_names);
prnt( '\n number of darks found : %g', num_dark)

% Projection range to read
if isempty( proj_range )
    proj_range = 1;
end
if numel( proj_range ) == 1
    proj_range = 1:proj_range:num_proj_found;
end
num_proj_used = numel( proj_range );
proj_nums = CellString2Vec( proj_names(proj_range), imtype_str_flag);
prnt( '\n number of projections found : %g', num_proj_found)
prnt( '\n number of projections used : %g', num_proj_used)
prnt( '\n projection range used : first:stride:last =  %g:%g:%g', proj_range(1), proj_range(2) - proj_range(1), proj_range(end))

%%% Read log files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new hdf5 log from statussever
h5log = sprintf('%s%s_nexus.h5', scan_path, scan_name);
im_shape_raw = [];
dtype = '';
tif_info = [];
offset_shift = 0;
% old scan-log, still needed
str = dir( sprintf( '%s*scan.log', scan_path) );
filename = sprintf( '%s/%s', str.folder, str.name);
[par, cur, cam] = p05_log( filename );
if isempty( eff_pixel_size )
    eff_pixel_size = par.eff_pixel_size;
end
eff_pixel_size_binned = raw_bin * eff_pixel_size;
if isempty( energy )
    if isfield( par, 'energy')
        energy = par.energy;
    end
end
if isempty( sample_detector_distance )
    sample_detector_distance = par.sample_detector_distance;
end
if ~exist( h5log, 'file')
    % Image shape and ROI
    filename = sprintf('%s%s', scan_path, ref_names{1});
    if ~raw_data
        [im_raw, tif_info] = read_image( filename );
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
        im_raw = read_raw( filename, im_shape_raw, dtype );
    end        
else
    % HDF5 log
    h5i = h5info( h5log );
    % energy, exposure time, image shape
    switch lower( cam )
        %%% CHECK h5 entry of camera1 / camera2 !!!!!!!!!!!!!!!!!!!!!!!!
        case 'ehd'
            energy = double(h5read( h5log, '/entry/hardware/camera1/calibration/energy') );
            exp_time = double(h5read( h5log, '/entry/hardware/camera1/calibration/exptime') );
            im_shape_raw = [3056 3056];
            dtype = 'uint16';
        case 'kit'
            energy = double(h5read( h5log, '/entry/hardware/camera2/calibration/energy') );
            exp_time = double(h5read( h5log, '/entry/hardware/camera2/calibration/exptime') );
            im_shape_raw = [5120 3840];
            dtype = 'uint16';
    end
    % Image shape
    filename = sprintf('%s%s', scan_path, ref_names{1});
    [im_raw, tif_info] = read_image( filename, '', [], tif_info, im_shape_raw, dtype );    
    % images
    stimg_name.value = unique( h5read( h5log, '/entry/scan/data/image_file/value') );
    stimg_name.time = h5read( h5log,'/entry/scan/data/image_file/time');
    stimg_key.value = h5read( h5log,'/entry/scan/data/image_key/value');
    stimg_key.time = double( h5read( h5log,'/entry/scan/data/image_key/time') );
    % PETRA ring current
    [petra.time, index] = unique( h5read( h5log,'/entry/hardware/beam_current/current/time') );
    petra.current = h5read( h5log,'/entry/hardware/beam_current/current/value');
    petra.current = petra.current(index);
    % rotation axis and wiggle
    s_rot.time = double( h5read( h5log, '/entry/scan/data/s_rot/time') );
    s_rot.value = h5read( h5log, '/entry/scan/data/s_rot/value');    
    s_stage_x.time = double( h5read( h5log, '/entry/scan/data/s_stage_x/time') );
    s_stage_x.value = h5read( h5log, '/entry/scan/data/s_stage_x/value');
    % wiggle di wiggle
    offset_shift = s_stage_x.value( ~boolean( stimg_key.value(par.n_dark+1:end) ) ) * 1e-3 / eff_pixel_size_binned;
    offset_shift = offset_shift - mean( offset_shift );
    offset_shift = offset_shift(proj_range);
    if visual_output(1) && numel( offset_shift ) > 2
        figure('Name', 'Rotation axis offset shift');
        plot( offset_shift, '.')
        title(sprintf('Rotation axis offset shift (zero mean)') )
        axis equal tight
        drawnow
    end
    % ring current
    X = double( petra.time(2:end) ); % first value is zero
    V = double( petra.current(2:end) ); % first value is zero
    Xq = double( stimg_name.time );
    stimg_name.current = (interp1( X, V, Xq, 'next', 100) + interp1( X, V, Xq + exp_time, 'previous', 100) ) / 2;
    cur_ref_val = stimg_name.current( stimg_key.value == 1 );
    cur_ref_name = stimg_name.value( stimg_key.value == 1 );    
    for nn = numel( cur_ref_name ):-1:1
        cur.ref(nn).val = cur_ref_val(nn);
        cur.ref(nn).name = cur_ref_name{nn};
        switch imtype_str_flag
            case 0
                cur.ref(nn).ind = str2double(cur.ref(nn).name(end-12:end-8));
            case 1
                cur.ref(nn).ind = str2double(cur.ref(nn).name(end-7:end-4));
        end
    end
    cur_proj_val = stimg_name.current( stimg_key.value == 0);
    cur_proj_name = stimg_name.value( stimg_key.value == 0);
    for nn = numel( cur_proj_name ):-1:1
        cur.proj(nn).val = cur_proj_val(nn);
        cur.proj(nn).name = cur_proj_name{nn};
        switch imtype_str_flag
            case 0                
                cur.proj(nn).ind = str2double(cur.proj(nn).name(end-12:end-8));
            case 1
                cur.proj(nn).ind = str2double(cur.proj(nn).name(end-7:end-4));
        end
    end           
    tomo.rot_angle.full_range = pi; % FIX/CHECK?
end
if ~isempty( raw_roi ) && ~isscalar( raw_roi ) && raw_roi(2) < 1
    raw_roi(2) = im_shape_raw(2) + raw_roi(2);
elseif isscalar(raw_roi) && raw_roi(1) < 1
    % Make roi_fac dependent on dark field
    if raw_roi(1) == -1
        roi_fac = 2;
    elseif raw_roi(1) < -1        
        roi_fac = abs( raw_roi );
    end
    % Read first non-zero flat
    mm = 1;
    while mean2( im_raw ) == 0
        mm = mm + 1;
        filename = sprintf('%s%s', scan_path, ref_names{mm});
        im_raw = read_image( filename, '', [], tif_info, im_shape_raw, dtype );    
    end
    % Read last non-zero flat
    mm = num_ref_found;
    filename = sprintf('%s%s', scan_path, ref_names{mm});
    im_raw2 = read_image( filename, '', [], tif_info, im_shape_raw, dtype );
    while mean2( im_raw2 ) == 0
        mm = mm - 1;
        filename = sprintf('%s%s', scan_path, ref_names{mm});
        im_raw2 = read_image( filename, '', [], tif_info, im_shape_raw, dtype );
    end
    % Compute crop position
    mm = median( (FilterPixel( im_raw, [0.1 0.1]) + FilterPixel( im_raw2, [0.1 0.1]) ) / 2, 1);    
    roi_thresh = roi_fac * double( min(mm) );
    [~,p] = min(abs(mm - roi_thresh));
    %%% IMPROVE STABILITY !!!
    dp = 200;
    if p < im_shape_raw(1)/2
        pl = p;
        [~,pr] = min(abs(mm(p+dp:end) - roi_thresh));
        pr = p + dp + pr;
    else
        pr = p;
        [~,pl] = min(abs(mm(1:p-dp) - roi_thresh));        
    end
    raw_roi = [pl pr];
    if visual_output(1)
        figure('Name', 'Auto ROI: raw image and cropping region');
        
        subplot(1,2,1)
        im = (single(im_raw) + single(im_raw2) ) / 2;
        pll = max( 1, pl - 20 );
        prr = min( im_shape_raw(2), pr + 20);
        im(:,[pll:pl, pr:prr] ) = max( im(:) );        
        %imsc1( im );
        imsc(im);
        title(sprintf('raw flat field: first + last'))        
        axis equal tight
        xticks('auto'),yticks('auto')
        camroll(90)
        
        subplot(1,2,2)
        plot( [mm'  roi_thresh*(ones(numel(mm),1))])
        camroll(90)
        title(sprintf('horizontal projection and cut level'))
        axis tight
        ax = gca;
        ax.YAxisLocation = 'right';
        ax.XDir = 'reverse';
        ax.YDir = 'normal';
        text(raw_roi(1), 2 * roi_thresh,sprintf('raw roi(1)=%u', raw_roi(1)))
        text(raw_roi(2), 2 * roi_thresh,sprintf('raw roi(2)=%u', raw_roi(2)))
        
        drawnow
    end
end
im_roi = read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype );
im_shape_binned = floor( size( im_roi) / raw_bin );
im_shape_binned1 = im_shape_binned(1);
im_shape_binned2 = im_shape_binned(2);
prnt( '\n energy : %.1f keV', energy / 1e3 )
prnt( '\n distance sample dector : %.1f mm', sample_detector_distance * 1000 )
prnt( '\n effective pixel size unbinned : %.2f micron',  eff_pixel_size * 1e6)
prnt( '\n raw image shape : %u x %u = %u pixels', im_shape_raw, prod( im_shape_raw ))
prnt( '\n raw image shape roi : %u x %u = %u pixels', size( im_roi ), numel( im_roi ) )
numel_im_roi_binned = prod( im_shape_binned );
prnt( '\n raw image shape roi binned : %u x %u = %u pixels', im_shape_binned, numel_im_roi_binned)
prnt( '\n raw binning factor : %u', raw_bin)
prnt( '\n raw binning before pixel filtering : %u', bin_before_filtering)

%%% Start parallel CPU pool %%%
t = toc;
[poolobj, poolsize] = OpenParpool(poolsize, use_cluster, [beamtime_path filesep 'scratch_cc']);
PrintVerbose( verbose && (poolsize > 1), '\nParpool opened on %s using %u workers. ', poolobj.Cluster.Profile, poolobj.Cluster.NumWorkers)
PrintVerbose( verbose && (poolsize > 1), ' done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 );
cdscandir = cd( scan_path );

%%% Read flat corrected projection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if read_flatcor(1)
    t = toc;
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
    prnt( '\n Read flat corrected projections.')
    if num_proj_read ~= num_proj_used
        fprintf('\n WARNING: Number of flat corrected projections read (%g) differs from number of projections to be processed (%g)!\n', num_proj_read, num_proj_used)
    end
    
    % Preallocation
    proj = zeros( im_shape_binned(1), im_shape_binned(2), num_proj_read, 'single');
    prnt( ' Allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
    proj_names_mat = NameCellToMat( proj_names );
    
    % Read flat corrected projections
    parfor nn = 1:num_proj_read
        filename = sprintf('%s%s', flatcor_path, proj_names_mat(nn, :));
        proj(:, :, nn) = fliplr( read_image( filename ) );
    end
    prnt( ' done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    
%%% Read and preprocess raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~read_flatcor    
    %%% Dark field %%
    t = toc;
    prnt( '\nProcessing dark fields.')
    darks = zeros( [im_shape_binned, num_dark], 'single');
    prnt( ' Allocated memory: %.2f MiB,', Bytes( darks, 2 ) )
    parfor nn = 1:num_dark
        filename = sprintf('%s%s', scan_path, dark_names{nn});
        im = single( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype) );
        % Remove large outliers. Assume Poisson distribtion at large lambda
        % is approximately a Gaussian distribution and set all value above
        % mean + 4 * std (99.994 of values lie within 4 std). Due to
        % outliers 4*std will contain much more values and is a good
        % estimate
        im_mean = mean( im(:) );
        im_std = std( im(:) );
        im( im > im_mean + 4*im_std) = im_mean;
        if bin_before_filtering
            darks(:, :, nn) = FilterPixel( Binning(im, raw_bin) / raw_bin^2, pixel_filter_threshold_dark );
        else
            darks(:, :, nn) = Binning( FilterPixel( im, pixel_filter_threshold_dark), raw_bin) / raw_bin^2;
        end
    end
    darks_min = min( darks(:) );
    darks_max = max( darks(:) );
    % Reject dark images which are all zero
    darks_to_use = zeros( 1, num_dark, 'logical' );
    parfor nn = 1:num_dark
        darks_to_use(nn) = boolean( max2( darks(:,:,nn) )  );
    end
    % Median dark
    dark = squeeze( median(darks(:,:,darks_to_use), 3) );
    %clear darks
    dark_med_min = min( dark(:) );
    dark_med_max = max( dark(:) );
    prnt( ' done in %.1f s', toc-t)
    if visual_output(1)
        h1 = figure('Name', 'data and flat-and-dark-field correction');
        subplot(2,3,1)
        imsc1( dark );        
        title(sprintf('median dark field'))
        colorbar
        axis equal tight
        xticks('auto'),yticks('auto')        
        drawnow
    end
    
    %%% Flat field %%%
    t = toc;
    prnt( '\nProcessing flat fields.')
    % Preallocation
    flat = zeros( [im_shape_binned, num_ref_used], 'single');
    num_zeros = zeros( 1, num_ref_used );
    prnt( ' Allocated memory: %.2f GiB,', Bytes( flat, 3 ) )
    % Parallel loop
    refs_to_use = zeros( 1, size( flat,3), 'logical');
    parfor nn = 1:num_ref_used
        filename = sprintf('%s%s', scan_path, ref_names_mat(nn,:));
        if bin_before_filtering
            flat(:, :, nn) = FilterPixel( Binning( single( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype ) ), raw_bin) / raw_bin^2, pixel_filter_threshold_flat);
        else
            flat(:, :, nn) = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype ), pixel_filter_threshold_flat), raw_bin) / raw_bin^2;
        end        
        % Check for zeros
        num_zeros(nn) =  sum( sum( flat(:,:,nn) < 1 ) );
        % Discard if any pixel is zero. Or should all pixels be zeros?
        refs_to_use(nn) = ~boolean( num_zeros(nn)  );        
    end
    
    % Delete empty refs
    flat(:,:,~refs_to_use) = [];
    
    % min/max values before dark field subtraction and ring current normalization
    flat_min = min( flat(:) );
    flat_max = max( flat(:) );
    
    % Dark field correction
    flat = bsxfun( @minus, flat, dark );
    
    % Ring current normalization
    if ring_current_normalization(1)                                            
        ref_ind_from_filenames = ref_nums;
        ref_ind_from_log = [cur.ref(ref_range).ind];
        if isequal( ref_ind_from_filenames, ref_ind_from_log )
            ref_rc = [cur.ref(ref_range).val];
            ref_rcm = mean( ref_rc(:) );
            scale_factor = 100 ./ shiftdim( ref_rc(refs_to_use), -1 );
            flat = bsxfun( @times, flat, scale_factor );
            if visual_output(1)
                hrc = figure('Name', 'Ring currents');
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
    
    % Enforce positivity
    parfor nn = 1:size( flat, 3 )
        im = flat(:,:,nn);
        m = im < 1;
        im(m) = 1;
        flat(:,:,nn) = im;
    end
    flat_min2 = min( flat(:) );
    flat_max2 = max( flat(:) );
    nn =  sum( flat(:) < 1 );
    if nn > 0
        fprintf('\n WARNING: flat field contains %u zeros\n', nn)
    end
    
    nn = sum( ~refs_to_use(:) );
    num_ref_used = num_ref_used - nn;
    prnt( ' done in %.1f s', toc-t)
    PrintVerbose( verbose && nn,'\n discarded empty refs : %u, %.2f%%', nn, 100*nn/num_ref_found )
    if sum( num_zeros )
        prnt( '\n flat fields with zeros :' )
        % print #zeros if not all pixels are zero
        for nn = 1:numel(num_zeros)
            if num_zeros(nn) ~= 0
                if isequal( num_zeros(nn), numel_im_roi_binned )
                    prnt( ' %u', nn )
                else
                    prnt( ' %u:%u', nn, num_zeros(nn) )
                end
            end
        end
    end
    
    % Show flat field
    if visual_output(1)
        if exist( 'h1' , 'var' )
            figure(h1)
        else
            h1 = figure('Name', 'data and flat-and-dark-field correction');
        end
        subplot(2,3,2)
        imsc1( flat(:,:,1) )        
        title(sprintf('flat field #1'))
        colorbar
        axis equal tight
        drawnow
    end
    
    %%% Projections %%%
    t = toc;
    prnt( '\nProcessing projections. ')
    % Preallocation
    proj = zeros( im_shape_binned(1), im_shape_binned(2), num_proj_used, 'single');
    prnt( ' Allocated memory: %.2f GiB,', Bytes( proj, 3 ) )
    img_names_mat = NameCellToMat( proj_names(proj_range) );
    
    % Display first raw image
    if visual_output(1)
        if exist( 'h1' , 'var' )
            figure(h1)
        else
            h1 = figure('Name', 'data and flat-and-dark-field correction');
        end
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        raw1 = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype ), pixel_filter_threshold_proj), raw_bin) / raw_bin^2;
        subplot(2,3,3)
        imsc1( raw1 )        
        title(sprintf('raw proj #1'))
        colorbar
        axis equal tight
        drawnow        
    end
    
    % Get absolut filter thresholds from percentage-wise pixel filtering
    % of 1st, middle, and last projection to speed up processing
    if pixel_filter_threshold_proj(1) < 1 || pixel_filter_threshold_proj(2) < 0.5        
        filename = sprintf('%s%s', scan_path, img_names_mat(num_proj_used, :));
        [~, ht(3), dt(3)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype ), pixel_filter_threshold_proj);        
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        [~, ht(2), dt(2)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype ), pixel_filter_threshold_proj);        
        filename = sprintf('%s%s', scan_path, img_names_mat(round(num_proj_used/2), :));
        [~, ht(1), dt(1)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype ), pixel_filter_threshold_proj);        
        HotThresh = median( ht );
        DarkThresh = median( dt );
    else
        HotThresh = pixel_filter_threshold_proj(1);
        DarkThresh = pixel_filter_threshold_proj(2);
    end
    
    % Read raw projections
    projs_to_use = zeros( 1, size( proj,3), 'logical' );
    num_zeros = zeros( 1, num_proj_used );
    parfor nn = 1:num_proj_used
        % Read and filter raw projection
        filename = sprintf('%s%s', scan_path, img_names_mat(nn,:));
        if bin_before_filtering
            im = FilterPixel( Binning( single( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype ) ), raw_bin) / raw_bin^2, [HotThresh, DarkThresh]);
        else
            im = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info, im_shape_raw, dtype ), [HotThresh, DarkThresh]), raw_bin) / raw_bin^2;
        end        
        % Check for zeros and reject images which are all zero
        num_zeros(nn) =  sum( sum( im < 1 ) );
        projs_to_use(nn) = ~boolean( num_zeros(nn)  );
        if projs_to_use(nn)
            proj(:, :, nn) = im;
        end
    end
    
    % Delete empty projections
    proj(:,:,~projs_to_use) = [];
    if offset_shift ~= 0
        offset_shift(~projs_to_use) = [];
    end
    
    raw_min = min( proj(:) );
    raw_max = max( proj(:) );    
    
    % Dark field correction
    proj = bsxfun( @minus, proj, dark);
    
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
            if visual_output(1)
                if exist( 'hrc', 'var' )
                    figure(hrc)
                else
                    hrc = figure('Name', 'Ring currents');
                end
%                 subplot(2,1,2);
%                 plot( proj_rc(:), '.' )
%                 axis tight
%                 title(sprintf('ring current: raw projections'))
%                 legend( sprintf( 'mean: %.2f mA', proj_rcm) )
%                 drawnow
                subplot(1,1,1);
                plot( ref_nums, ref_rc(:), '.',proj_nums, proj_rc(:), '.' )
                axis tight
                title(sprintf('ring currents'))
                legend( sprintf( 'flats, mean: %.2f mA', ref_rcm), sprintf( 'projs, mean: %.2f mA', proj_rcm) )
                drawnow
            end
        else
            fprintf('\n WARNING: projections not normalized by ring current. Names read from directory and log-file are inconsistent.\n')
        end
    end
    
    % Enforce positivity
    parfor nn = 1:size( proj, 3 )
        im = proj(:,:,nn);
        m = im < 1;
        if sum( m(:) ) > 0
            im(m) = 1;
            proj(:,:,nn) = im;
        end
    end
    
    raw_min2 = min( proj(:) );
    raw_max2 = max( proj(:) );
    num_empty = sum( ~projs_to_use(:) );
    num_proj_used = num_proj_used - num_empty;
    prnt( ' done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    PrintVerbose( verbose && num_empty,'\n discarded empty projections : %u, %.2f%%', num_empty, 100*num_empty/size(proj,3) )
    if sum( num_zeros )
        prnt( '\n projections with zeros :' )
        % print #zeros if not all pixels are zero
        for nn = 1:numel(num_zeros)
            if num_zeros(nn) ~= 0
                if isequal( num_zeros(nn), numel_im_roi_binned )
                    prnt( ' %u', nn )
                else
                    prnt( ' %u:%u', nn, num_zeros(nn) )
                end
            end
        end
    end
    prnt( '\n hot- / dark-pixel filter threshold : %f, %f', HotThresh, DarkThresh )
    prnt( '\n global min/max of raws after filtering and binning:  %6g %6g', raw_min, raw_max)
    prnt( '\n global min/max of raws after dark-field correction and ring current normalization:  %6g %6g', raw_min2, raw_max2)
    
    if stop_after_data_reading
        fprintf( '\n' )
        keyboard;
    end
    
    %% Projection/flat field correlation and flat field correction %%%%%%%%
    %%%% STOP HERE FOR FLATFIELD CORRELATION MAPPING %%%%%%%%%%%%%%%%%%%%%%
    [proj, corr, proj_roi, flat_roi] = proj_flat_correlation( proj, flat, image_correlation.method, image_correlation.area_width, image_correlation.area_height, im_shape_binned, image_correlation.shift.max_pixelshift, image_correlation.num_flats, decimal_round_precision, flatcor_path, verbose, visual_output);
    proj_min0 = min( proj(:) );
    proj_max0 = max( proj(:) );
    prnt( '\n global min/max after flat-field corrected:  %6g %6g', proj_min0, proj_max0);   

    % Experimental: Filter strong absorption
    if strong_abs_thresh < 1
        t = toc;
        prnt( '\n set flat-corrected values below %f to one, ', strong_abs_thresh)
        parfor nn = 1:size( proj, 3 )
            im = proj(:,:,nn);
            m = im < strong_abs_thresh;
            im(m) = 1;
            proj(:,:,nn) = im;
        end
            prnt( ' done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    end
    
    prnt( '\n sinogram size = [%g, %g, %g]', size( proj ) )
    if visual_output(1)
        if exist( 'h1' , 'var' )
            figure(h1)
        else
            h1 = figure('Name', 'data and flat-and-dark-field correction');
        end
        
        subplot(2,3,4)
        imsc1( proj(:,:,1))        
        xticks([]);yticks([])
        title(sprintf('intensity: first proj'))
        colorbar
        axis equal tight
        
        subplot(2,3,5)
        imsc1( proj(:,:,round(size(proj,3)/2)))        
        xticks([]);yticks([])
        title(sprintf('intensity: middle proj'))
        colorbar
        axis equal tight
        
        subplot(2,3,6)
        imsc1( proj(:,:,end))        
        xticks([]);yticks([])
        title(sprintf('intensity: last proj'))
        colorbar
        axis equal tight
        
        drawnow
    end
    
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
        angles = s_rot.value( ~boolean( stimg_key.value(par.n_dark+1:end) ) ) * pi / 180;        
        angles = angles(proj_range);
    else        
        num_proj = par.num_proj;
        switch lower( cam )
            case 'ehd'
                angles = tomo.rot_angle.full_range * (0:num_proj - 1) / (num_proj - 1); % EHD: ok
            case 'kit'
                angles = tomo.rot_angle.full_range * (0:num_proj - 1) / num_proj; % KIT: ok if par.projections exist
        end
    end
    % drop angles where projections are empty
    angles(~projs_to_use) = [];

    %% Ring artifact filter    
    if ring_filter.apply && ring_filter.apply_before_stitching
        [proj, h] = p05_filter_ring_artefacts( ring_filter, proj, angles, verbose, visual_output);
    end
    proj_min = min( proj(:) );
    proj_max = max( proj(:) );
    
    %% Write corrected projections
    if write.flatcor
        t = toc;
        prnt( '\nSave flat-corrected projections.')
        CheckAndMakePath( flatcor_path )
        % Projections
        parfor nn = 1:size( proj, 3 )
            filename = sprintf('%sproj_%06u.tif', flatcor_path, nn );
            write32bitTIFfromSingle(filename, rot90( proj(:, :, nn) ) );
        end
        prnt( ' done in %.1f (%.2f min)', toc-t, (toc-t)/60)
    end
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if phase_retrieval.apply
    if isempty( tomo.take_neg_log )
        tomo.take_neg_log = 0;
    end
    if phase_retrieval.apply_before
        edp = [energy, sample_detector_distance, eff_pixel_size_binned];
        [proj, reco_phase_path] = p05_phase_retrieval( phase_retrieval, edp, proj, write.path, write.subfolder.reco, write, phase_map_path, verbose, visual_output);
    end
end

%%% Rotation axis position and tomgraphic reconstruction parameters %%%
tint = 0;
if tomo.run || tomo.run_interactive_mode
    prnt( '\nTomography:')
    t = toc;
    % ROI for correlation of projections at angles 0 & pi
    if isempty( tomo.rot_axis.corr_area1 )
        switch excentric_rot_axis
            case -1
                tomo.rot_axis.corr_area1 = [0 0.25];
            case 0
                tomo.rot_axis.corr_area1 = [0.25 0.75];
            case 1
                tomo.rot_axis.corr_area1 = [0.75 1];
        end
    end
    if isempty( tomo.rot_axis.corr_area2 )
        tomo.rot_axis.corr_area2 = [0.1 0.9];
    end
    tomo.rot_axis.corr_area1 = IndexParameterToRange( tomo.rot_axis.corr_area1, im_shape_binned1 );
    tomo.rot_axis.corr_area2 = IndexParameterToRange( tomo.rot_axis.corr_area2, im_shape_binned2 );
    
    % Full rotation angle
    if isempty( tomo.rot_angle.full_range )
        if isfield( par, 'rotation')
            % From log file
            tomo.rot_angle.full_range = par.rotation / 180 * pi;
        elseif exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
            % from beam current log
            tomo.rot_angle.full_range = (cur.proj(end).angle - cur.proj(1).angle) * pi /180; % KIT: , EHD: ok
        else
            % Guess from correlation of projections
            im1 = proj( tomo.rot_axis.corr_area1, tomo.rot_axis.corr_area2, 1);
            im2 = proj( tomo.rot_axis.corr_area1, tomo.rot_axis.corr_area2, end);
            [~, cm1] = ImageCorrelation( im1, im2, 0, 0, 0); % large if full angle is  2 * pi
            [~, cm2] = ImageCorrelation( im1, flipud(im2), 0, 0, 0); % large if full angle is pi
            if max( cm1(:) ) > max( cm2(:) )
                tomo.rot_angle.full_range = 2 * pi;
            elseif max( cm1(:) ) < max( cm2(:) )
                tomo.rot_angle.full_range = pi;
            else
                error('Determination of full angle of rotation from correlation of first and last projection not successful.')
            end
        end
    end
    prnt( '\n full rotation angle: %g * pi', tomo.rot_angle.full_range / pi)
    if numel( angles ) ~= num_proj_used
        error('Number of elements in array of angles (%g) unequal number of projections read (%g)', numel( angles ), num_proj_used)
    end
    % retrieve index at angles 0 and pi
    [val1, ind1] = min( abs( angles ));
    [val2, ind2] = min( abs( angles - pi ));
    % ROI
    im1 = proj( tomo.rot_axis.corr_area1, tomo.rot_axis.corr_area2, ind1);
    im2 = flipud( proj( tomo.rot_axis.corr_area1, tomo.rot_axis.corr_area2, ind2) );
    prnt( '\n correlated images : [filename index, projection index, angle / pi] = [%g, %g, %g] and [%g, %g, %g]', ...
        proj_nums(ind1), ind1, val1 / pi, proj_nums(ind2), ind2, (val2 + pi) / pi )
    if tomo.rot_axis.corr_gradient(1)
        l = 2;
        [g11, g12] = gradient(im1);
        im1 = abs(g11).^l + abs(g12).^l;
        [g21, g22] = gradient(im2);
        im2 = abs(g21).^l + abs(g22).^l;
    end
    % Correlation
    out = ImageCorrelation( im1, im2, 0, 0, 0);
    % relative shift
    rot_corr_shift_x = round( out.shift1, decimal_round_precision) + tomo.rot_axis.corr_area1(1) + (tomo.rot_axis.corr_area1(end) - im_shape_binned1) - 1;
    rot_corr_shift_y = round( out.shift2, decimal_round_precision) + tomo.rot_axis.corr_area2(1) + (tomo.rot_axis.corr_area1(end) - im_shape_binned1) - 1;
    prnt( '\n relative shift: %g, %g', rot_corr_shift_x, rot_corr_shift_y)
    rot_axis_offset_calc = rot_corr_shift_x / 2;
    if numel( offset_shift ) > 1
        rot_axis_offset_calc = offset_shift(ind1) - offset_shift(ind2);
    end
    rot_axis_pos_calc = im_shape_binned1 / 2 + rot_axis_offset_calc;
    prnt( '\n calculated rotation axis offset: %.2f', rot_axis_offset_calc)
    % use calculated offset if both offset and position are empty
    if isempty( tomo.rot_axis.offset ) && isempty( tomo.rot_axis.position )
        tomo.rot_axis.offset = rot_axis_offset_calc;
    end
    if isempty( tomo.rot_axis.position )
        tomo.rot_axis.position = im_shape_binned1 / 2 + tomo.rot_axis.offset;
    end
    if isempty( tomo.rot_axis.offset )
        tomo.rot_axis.offset = tomo.rot_axis.position - im_shape_binned1 / 2;
    end
    
    % Tilt of rotation axis
    if interactive_mode.rot_axis_tilt(1)
        %im1c = RotAxisSymmetricCropping( proj(:,tomo.rot_axis.corr_area2,ind1), tomo.rot_axis.position, 1);
        %im2c = flipud(RotAxisSymmetricCropping( proj(:,tomo.rot_axis.corr_area2,ind2) , tomo.rot_axis.position, 1));
        corr_offset = ( offset_shift(ind1) + offset_shift(ind2) ) / 2;
        im1c = RotAxisSymmetricCropping( proj(:,tomo.rot_axis.corr_area2,ind1), tomo.rot_axis.position + corr_offset, 1);
        im2c = flipud(RotAxisSymmetricCropping( proj(:,tomo.rot_axis.corr_area2,ind2) , tomo.rot_axis.position + corr_offset, 1));
        
        [optimizer, metric] = imregconfig('monomodal');
        tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);
        rot_axis_tilt_calc = asin( tform_calc.T(1,2) ) / 2;
        prnt( '\n calculated tilt of rotation axis : %g rad = %g deg', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
    else
        rot_axis_tilt_calc = [];
    end
    if isempty( tomo.rot_axis.tilt )
        tomo.rot_axis.tilt = 0;
    end    
        
    %% Interactive mode:rotation axis position / tilt %%%%%%%%%%%%%%%%%%%%%
    
%    %%% AUTOMATIC MODE %%%
%     automatic_mode = 0; % Find rotation axis position automatically. NOT IMPLEMENTED!
%     automatic_mode_coarse = 'entropy'; % NOT IMPLEMENTED!
%     automatic_mode_fine = 'iso-grad'; % NOT IMPLEMENTED!
%     if automatic_mode(1) % NOT IMPLEMENTED YET
%         pts = 10;
%         im_center = im_shape_binned1 / 2;
%         offset_stride = floor( im_center / pts );
%         offset = -im_center:offset_stride:im_center;        
%     end
    
    tint = 0;
    if interactive_mode.rot_axis_pos(1)
        tint = toc;
        fprintf( '\n\nENTER INTERACTIVE MODE' )
        fprintf( '\n number of pixels: %u', im_shape_binned1)
        fprintf( '\n image center: %.1f', im_shape_binned1 / 2)
        
        if phase_retrieval.apply && phase_retrieval.apply_before
            itake_neg_log = 0;
        else
            itake_neg_log = 1;
        end
        if phase_retrieval.apply
            inumber_of_stds = 9;
        else
            inumber_of_stds = 4;
        end
        if interactive_mode.slice_number > 1
            slice = interactive_mode.slice_number;
        elseif interactive_mode.slice_number <= 1 && interactive_mode.slice_number >= 0
            slice = round((im_shape_binned2 - 1) * interactive_mode.slice_number + 1 );
        end
        fprintf( '\n slice : %u', slice)        
        fprintf( '\n\nOFFSET:' )
        fprintf( '\n current rotation axis offset / position : %.2f, %.2f', tomo.rot_axis.offset, tomo.rot_axis.position)
        fprintf( '\n calcul. rotation axis offset / position : %.2f, %.2f', rot_axis_offset_calc, rot_axis_pos_calc)
        fprintf( '\n default offset range : current ROT_AXIS_OFFSET + (-4:0.5:4)')
        offset = input( '\n\nENTER RANGE OF ROTATION AXIS OFFSETS\n (if empty: use default range, scalar: skip interactive mode, ''s'': change slice, ''d'': debug mode): ');        
        if ~isempty( offset ) && strcmp( offset(1), 's' )
            slice = input( sprintf( '\n\nENTER ABSOLUTE [1,%u] OR RELATIVE [0,1] SLICE NUMBER : ', im_shape_binned2) );            
            if slice <= 1 && slice >= 0
                slice = round((im_shape_binned2 - 1) * slice + 1 );
            end
            fprintf( ' \n new slice : %u', slice );
            offset = input( '\n\nENTER RANGE OF ROTATION AXIS OFFSETS\n (if empty use default range, if scalar skips interactive mode): ');
        end
        if ~isempty( offset ) && strcmp( offset(1), 'd')
            keyboard
            offset = input( '\n\nENTER RANGE OF ROTATION AXIS OFFSETS\n (if empty use default range, if scalar skips interactive mode): ');
        end
        if isempty( offset )
            % default range is centered at the given or calculated offset
            offset = tomo.rot_axis.offset + (-4:0.5:4);
        end
        if isscalar( offset )
            fprintf( ' old rotation axis offset : %.2f', tomo.rot_axis.offset)
            tomo.rot_axis.offset = offset;
            fprintf( '\n new rotation axis offset : %.2f', tomo.rot_axis.offset)
        end
        
        % Loop over offsets
        while ~isscalar( offset )
            
            [ivol_shape, ivol_size] = volshape_volsize( proj, tomo.vol_shape, tomo.vol_size, median(offset), verbose);
            
            % Reco
            %%% Clean up !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            fra = tomo;
            fra.angles = angles;
            fra.slice = slice;
            fra.offset = offset;
            fra.tilt = tomo.rot_axis.tilt;
            fra.take_neg_log = itake_neg_log;
            fra.inumber_of_stds = inumber_of_stds;
            fra.vol_shape = ivol_shape;
            fra.vol_size = ivol_size;
            fra.lamino = interactive_mode.lamino;
            fra.fixed_tilt = interactive_mode.fixed_other_tilt;
            fra.astra_gpu_index = tomo.astra_gpu_index;
            fra.offset_shift = offset_shift;
            [vol, metrics_offset] = find_rot_axis_offset( fra, proj );
            
            % Metric minima
            [~, min_pos] = min(cell2mat({metrics_offset(:).val}));
            [~, max_pos] = max(cell2mat({metrics_offset(:).val}));
            
            % Print image number, rotation axis values, and different metrics
            fprintf( '\n\nOFFSET:' )
            fprintf( '\n current rotation axis offset/position : %.2f, %.2f', tomo.rot_axis.offset, tomo.rot_axis.position)
            fprintf( '\n current slice : %u\n', slice )
            fprintf( '%11s', 'image no.', 'offset', metrics_offset.name)
            for nn = 1:numel(offset)
                if offset(nn) == tomo.rot_axis.offset
                    cprintf( 'Green', sprintf('\n%11u%11.2f', nn, offset(nn)))
                else
                    cprintf( 'Black', '\n%11u%11.2f', nn, offset(nn))
                end
                
                for mm = 1:numel(metrics_offset)
                    if min_pos(mm) == nn
                        cprintf( 'Red', '%11.3g', metrics_offset(mm).val(nn) )
                    elseif max_pos(mm) == nn
                        cprintf( 'Blue', '%11.3g', metrics_offset(mm).val(nn) )
                    else
                        cprintf( 'Black', '%11.3g', metrics_offset(mm).val(nn) )
                    end
                end
            end
            
            % Plot metrics
            h_rot_off = figure('Name', 'OFFSET: metrics');
            x = 1:7;%[1:4 6:7];
            Y = cell2mat({metrics_offset(x).val});
            plot( offset, Y, '-+');            
            axis tight
            xlabel( 'offset' )
            legend( metrics_offset(x).name )
            ax1 = gca;            
            ax2 = axes( 'Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
            line(1:numel( offset ), 0, 'Parent', ax2 )
            xlabel( 'index (image no.)' )
            set( ax1, 'YTick', [] ) % 'XTickMode', 'auto', 'XMinorTick', 'on')
            set( ax2, 'YTick', [] )
            title(sprintf('rotation axis: metrics VS offset'))
            drawnow
            
            % Play
            nimplay(vol, 1, [], 'OFFSET: sequence of reconstructed slices using different rotation axis offsets')
            
            % Input
            offset = input( '\n\nENTER ROTATION AXIS OFFSET OR A RANGE OF OFFSETS\n (if empty use current offset, ''s'' to change slice number, ''d'' for debug mode): ');
            if ~isempty( offset ) && strcmp( offset(1), 's' )
                slice = input( sprintf( '\n\nENTER ABSOLUTE [1,%u] OR RELATIVE [0,1] SLICE NUMBER : ', im_shape_binned2) );
                if slice <= 1 && slice >= 0
                    slice = round((im_shape_binned2 - 1) * slice + 1 );
                end
                fprintf( ' \n new slice : %u', slice );
                offset = input( '\n\nENTER RANGE OF ROTATION AXIS OFFSETS\n (if empty use default range, if scalar skips interactive mode): ');
            end
            if ~isempty( offset ) && strcmp( offset(1), 'd')
                keyboard
                offset = input( '\n\nENTER RANGE OF ROTATION AXIS OFFSETS\n (if empty use default range, if scalar skips interactive mode): ');
            end
            if isempty( offset )
                offset = tomo.rot_axis.offset + (-4:0.5:4);
            end
            if isscalar( offset )
                fprintf( ' old rotation axis offset : %.2f', tomo.rot_axis.offset)
                tomo.rot_axis.offset = offset;
                fprintf( '\n new rotation axis offset : %.2f', tomo.rot_axis.offset)
                
                % TILT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if interactive_mode.rot_axis_tilt(1)
                    fprintf( '\n\nTILT:' )
                    fprintf( '\n current rotation axis tilt : %g rad = %g deg', tomo.rot_axis.tilt, tomo.rot_axis.tilt * 180 / pi)
                    fprintf( '\n calcul. rotation axis tilt : %g rad = %g deg', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
                    fprintf( '\n default tilt range is : current ROT_AXIS_TILT + (-0.005:0.001:0.005)')
                    tilt = input( '\n\nENTER TILT OF ROTATION AXIS OR RANGE OF TILTS\n (if empty: use default, ''s'': change slice, ''d'': debug mode):');
                    % option to change which slice to reconstruct
                    if ~isempty( tilt ) && strcmp( tilt(1), 's' )                                                
                        slice = input( sprintf( '\n\nENTER ABSOLUTE [1,%u] OR RELATIVE [0,1] SLICE NUMBER : ', im_shape_binned2) );
                        if slice <= 1 && slice >= 0
                            slice = round((im_shape_binned2 - 1) * slice + 1 );
                        end
                        tilt = input( '\n\nENTER TILT OF ROTATION AXIS OR RANGE OF TILTS\n (if empty use default):');
                    end
                    if ~isempty( tilt ) && strcmp( tilt(1), 'd')
                        keyboard;
                        tilt = input( '\n\nENTER TILT OF ROTATION AXIS OR RANGE OF TILTS\n (if empty use default):');
                    end
                    if isempty( tilt )
                        tilt = tomo.rot_axis.tilt + (-0.005:0.001:0.005);
                    end
                    
                    while ~isscalar( tilt )
                        
                        % Reco                        
                        fra.angles = angles;
                        fra.slice = slice;
                        fra.offset = offset;
                        fra.tilt = tilt;
                        fra.take_neg_log = itake_neg_log;
                        fra.inumber_of_stds = inumber_of_stds;
                        fra.vol_shape = ivol_shape;
                        fra.vol_size = ivol_size;
                        fra.lamino = interactive_mode.lamino;
                        fra.fixed_tilt = interactive_mode.fixed_other_tilt;
                        fra.astra_gpu_index = tomo.astra_gpu_index;
                        fra.offset_shift = offset_shift;
                        [vol, metrics_tilt] = find_rot_axis_tilt( fra, proj);
                        
                        % Metric minima
                        [~, min_pos] = min(cell2mat({metrics_tilt(:).val}));
                        [~, max_pos] = max(cell2mat({metrics_tilt(:).val}));
                        
                        % Print image number and rotation axis tilt
                        fprintf( '%11s', 'image no.', 'tilt/rad', 'tilt/deg', metrics_tilt.name )
                        for nn = 1:numel(tilt)
                            if tilt(nn) == tomo.rot_axis.tilt
                                cprintf( 'Green', sprintf( '\n%11u%11g%11g', nn, tilt(nn), tilt(nn)/pi*180 ) )
                            else
                                cprintf( 'Black', sprintf( '\n%11u%11g%11g', nn, tilt(nn), tilt(nn)/pi*180 ) )
                            end
                            for mm = 1:numel(metrics_tilt)
                                if min_pos(mm) == nn
                                    cprintf( 'Red', '%11.3g', metrics_tilt(mm).val(nn) )
                                elseif max_pos(mm) == nn
                                    cprintf( 'Blue', '%11.3g', metrics_tilt(mm).val(nn) )
                                else
                                    cprintf( 'Black', '%11.3g', metrics_tilt(mm).val(nn) )
                                end
                            end
                        end
                        
                        % Plot metrics
                        h_rot_tilt = figure('Name', 'TILT: metrics');
                        x = 6:7;
                        Y = cell2mat({metrics_tilt(x).val});
                        plot( tilt, Y, '-+');
                        axis tight
                        xlabel( 'tilt angle' )
                        legend( metrics_tilt(x).name)
                        ax1 = gca;
                        set( ax1, 'YTick', [] )
                        ax2 = axes( 'Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
                        line(1:numel( offset ), 0, 'Parent', ax2 )
                        xlabel( 'index (image no.)' )
                        set( ax2, 'YTick', [] )
                        title(sprintf('rotation axis: metrics VS tilt'))
                        drawnow
                        
                        % Play
                        nimplay(vol, 1, [], 'TILT: sequence of reconstructed slices using different rotation axis tilts')
                        
                        % Input
                        tilt = input( '\nENTER TILT OF ROTATION AXIS OR RANGE OF TILTS\n (if empty use current tilt):');
                        if isempty( tilt )
                            tilt = tomo.rot_axis.tilt;
                        end
                        
                        if isscalar( tilt )
                            tomo.rot_axis.tilt = tilt;
                            
                            tomo.rot_axis.position = im_shape_binned1 / 2 + tomo.rot_axis.offset;
                            
                            % Compare projection at 0 pi and projection at 1 pi corrected for rotation axis tilt
%                             im1c = RotAxisSymmetricCropping( proj(:,tomo.rot_axis.corr_area2,ind1), tomo.rot_axis.position, 1);
%                             im2c = flipud(RotAxisSymmetricCropping( proj(:,tomo.rot_axis.corr_area2,ind2) , tomo.rot_axis.position, 1));
                            corr_offset = ( offset_shift(ind1) + offset_shift(ind2) ) / 2;
                            im1c = RotAxisSymmetricCropping( proj(:,tomo.rot_axis.corr_area2,ind1), tomo.rot_axis.position + corr_offset, 1);
                            im2c = flipud(RotAxisSymmetricCropping( proj(:,tomo.rot_axis.corr_area2,ind2) , tomo.rot_axis.position + corr_offset, 1));
                            [optimizer, metric] = imregconfig('monomodal');
                            tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);                            
                            rot_axis_tilt_calc = asin( tform_calc.T(1,2) ) / 2;
                            im2c_warped_calc =  imwarp(im2c, tform_calc, 'OutputView', imref2d(size(im1c)));
                            
                            tform_int = tform_calc;
                            tform_int.T = [cos( 2 * tomo.rot_axis.tilt ) sin( 2 * tomo.rot_axis.tilt ) 0; ...
                                -sin( 2 * tomo.rot_axis.tilt ) cos( 2 * tomo.rot_axis.tilt ) 0 ; ...
                                tform_calc.T(3,1) tform_calc.T(3,2) 1];                            
                            % Translation: [1 0 0, 0 1 0, t_x t_y 0]
                            % Scale: [s_x 0 0, 0 s_y 0, 0 0 1]
                            % Shear: [1 sh_y 0, sh_x 1 0, 0 0 1]
                            % Rotation: [cos(t) sin(t) 0, -sin(t) cos(t) 0, 0 0 1]

                            im2c_warped_int =  imwarp(im2c, tform_int, 'OutputView', imref2d(size(im1c)));
                            
                            xt = ceil( 3 * abs( sin(2*tomo.rot_axis.tilt) ) * max( size(im1c)) ) + 2;
                            
                            fprintf( '\n current rotation axis tilt from interactive mode: %g rad (%g deg)', tomo.rot_axis.tilt, tomo.rot_axis.tilt * 180 / pi)
                            fprintf( '\n calcul. rotation axis tilt from registration    : %g rad (%g deg)', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
                            
                            name = sprintf( 'registered projections at %g and %g degree. rot axis tilt from INTERACTIVE mode: %g, rot axis offset: %g', angles(ind1)/pi*180, angles(ind2)/pi*180, tomo.rot_axis.tilt, tomo.rot_axis.offset);
                            nimplay( cat(3, im1c(xt:end-xt,xt:end-xt)', im2c_warped_int(xt:end-xt,xt:end-xt)'), 1, 0, name)
                            
                            name = sprintf( 'registered projections at %g and %g degree. corrected. rot axis tilt from REGISTRATION: %g, rot axis offset: %g', angles(ind1)/pi*180, angles(ind2)/pi*180, rot_axis_tilt_calc, tomo.rot_axis.offset);
                            nimplay( cat(3, im1c(xt:end-xt,xt:end-xt)', im2c_warped_calc(xt:end-xt,xt:end-xt)'), 1, 0, name)
                            
                            tilt = input( '\nENTER ROTATION AXIS TILT\n (if empty use current value): ');
                            if isempty( tilt )
                                tilt = tomo.rot_axis.tilt;
                            else
                                tomo.rot_axis.tilt = tilt;
                            end
                            
                            offset = input( '\nENTER RANGE OF OFFSETS TO CONTINUE INTERACTIVE LOOP OR TYPE ENTER TO EXIT LOOP: ');
                            if isempty( offset )
                                offset = tomo.rot_axis.offset;
                            end
                        end
                    end
                else
                    tomo.rot_axis.tilt = 0;
                end
            end
        end
        
        tomo.rot_axis.position = im_shape_binned1 / 2 + tomo.rot_axis.offset;
        fprintf( '\nEND OF INTERACTIVE MODE' )
        
        tint = toc - tint;
    end
    prnt( '\n rotation axis offset: %.2f', tomo.rot_axis.offset );
    prnt( '\n rotation axis position: %.2f', tomo.rot_axis.position );
    prnt( '\n rotation axis tilt: %g rad (%g deg)', tomo.rot_axis.tilt, tomo.rot_axis.tilt * 180 / pi)
    [tomo.vol_shape, tomo.vol_size] = volshape_volsize( proj, tomo.vol_shape, tomo.vol_size, tomo.rot_axis.offset, verbose);
    
    if interactive_mode.rot_axis_tilt(1) && visual_output(1)
        h4 = figure('Name','Projections at 0 and pi cropped symmetrically to rotation center');
        n = 2;
        m = 2;
        
        subplot(m, n, 1)
        imsc1( im1c )
        axis equal tight
        title(sprintf('proj at 0'))
        colorbar
        
        subplot(m, n, 2)
        imsc1( im2c )
        axis equal tight
        title(sprintf('proj at pi'))
        colorbar
                
        xt = ceil( 3 * abs( sin(2*tomo.rot_axis.tilt) ) * max( size(im1c)) ) + 2;
        subplot(m, n, 3)
        imsc1( abs( im1c(xt:end-xt,xt:end-xt) - im2c(xt:end-xt,xt:end-xt) ) )
        axis equal tight
        title(sprintf('difference before'))
        colorbar
        
        subplot(m, n, 4)
        
        %%% tform_cal or tform_int !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        [optimizer, metric] = imregconfig('monomodal');
        tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);
        im2c_warped_int =  imwarp(im2c, tform_calc, 'OutputView',imref2d(size(im1c)));
        imsc1( abs( im1c(xt:end-xt,xt:end-xt) - im2c_warped_int(xt:end-xt,xt:end-xt) ) )
        axis equal tight
        title(sprintf('difference corrected'))
        colorbar
        
        drawnow
    end
end

%% Stitch projections
if stitch_projections(1)
    t = toc;
    prnt( '\nStitch projections:')
    if tomo.rot_angle.full_range < 1.9 * pi
        error( 'full angle of rotation smaller than 2 pi: %g pi', tomo.rot_angle.full_range/pi)
    end
    % number of unstitched projections
    num_proj = size( proj, 3);
    % last projection within [0,pi)
    [~, num_proj_sti] = min( abs(angles - pi));
    % number of stitched projections
    num_proj_sti = num_proj_sti - 1;
    % index range of projections to be stitched
    xl = 1:round(tomo.rot_axis.position);
    xr = 1:xl(end)-1;
    im_shape_sti1 = numel( xl ) + numel( xr );
    % Preallocation
    proj_sti = zeros( im_shape_sti1 , im_shape_binned2, num_proj_sti, 'single');
    for nn = 1:num_proj_sti
        nn2 = mod(num_proj_sti + nn - 1, num_proj) + 1;
        im = zeros( im_shape_sti1, im_shape_binned2);
        switch lower( stitch_method )
            case 'step'
                im = cat(1, proj(xl,:,nn), flipud( proj(xr,:,nn2) ) );
            case {'linear', 'sine'}
                % overlap region
                overlap = round(2 * tomo.rot_axis.position) - im_shape_binned1 : im_shape_binned1;
                % overlap ramp
                x = (0:1/(numel(overlap)-1):1);
                % 1D weight
                w = ones(im_shape_binned1, 1);
                switch lower( stitch_method )
                    case 'linear'
                        w(overlap) = 1 - x;
                    case 'sine'
                        w(overlap) = 0.5 * cos(pi*x) + 0.5;
                end
                % weighted projections
                iml = bsxfun(@times, proj(:,:,nn), w);
                imr = flipud( bsxfun(@times, proj(:,:,nn2), w ) );
                % stitched projection
                im(1:im_shape_binned1,:) = iml;
                im(end - im_shape_binned1 + 1:end,:) = im(end - im_shape_binned1 + 1:end,:) + imr;
        end
        proj_sti(:,:,nn) = im;
    end
    pause(0.01)
    proj = proj_sti;
    clear proj_sti;
    angles = angles(1:num_proj_sti);
    prnt( ' done in %.1f (%.2f min)', toc-t, (toc-t)/60)
    prnt( '\n shape of stitched projections : %u %u %u', size( proj ) )
    prnt( '\n memory allocated : %.2f GiB', Bytes( proj, 3 ) )
end

%% Ring artifact filter %%
if ring_filter.apply && ~ring_filter.apply_before_stitching    
    [proj, h] = p05_filter_ring_artefacts( ring_filter, proj, angles, verbose, visual_output);
    proj_min = min( proj(:) );
    proj_max = max( proj(:) );
end

%% Crop projections at rotation axis position %%
if crop_at_rot_axis(1)
    % This is to avoid oversampling for scans with excentric rotation axis
    % and reconstructing WITHOUT stitching
    switch excentric_rot_axis
        case 1
            proj( ceil(tomo.rot_axis.position) + 1:end, :, :) = [];
        case -1
            %%% CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            proj( 1:floor(tomo.rot_axis.position)-1, :, :) = [];
    end
    if isempty( vol_shape )
        tomo.vol_shape = [raw_im_shape_binned1, raw_im_shape_binned1, raw_im_shape_binned2];
    end
end

%% Save sinograms %%
if write.sino
    t = toc;
    prnt( '\nSave sinogram:')
    CheckAndMakePath(sino_path)
    % Save slices
    parfor nn = 1:im_shape_binned2
        filename = sprintf( '%ssino_%06u.tif', sino_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    pause(0.01)
    prnt( ' done in %.1f s (%.2f min)', toc-t, (toc-t)/60)
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if phase_retrieval.apply && ~phase_retrieval.apply_before
    edp = [energy, sample_detector_distance, eff_pixel_size_binned];    
    [proj, reco_phase_path] = p05_phase_retrieval( phase_retrieval, edp, proj, write.path, write.subfolder.reco, write, phase_map_path, verbose, visual_output);
end

%%% Save sinograms of phase maps %%%
if phase_retrieval.apply && write.phase_sino
    t = toc;
    prnt( '\nSave phase map sinogram:')
    CheckAndMakePath(sino_phase_path)
    % Save slices
    parfor nn = 1:im_shape_binned2
        filename = sprintf( '%ssino_%06u.tif', sino_phase_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    pause(0.01)    
    prnt( ' done in %.1f s (%.2f min)', toc - t, (toc - t) / 60)
end

%%% Post phase retrieval binning %%%
phase_bin = phase_retrieval.post_binning_factor; % alias for readablity
if phase_retrieval.apply && ( phase_bin(1) > 1)
    t = toc;
    prnt( '\nPost phase retrieval binning:')
    proj_bin = zeros( floor(size( proj ) ./ [phase_bin phase_bin 1] ), 'single');
    parfor nn = 1:size( proj, 3)
        proj_bin(:,:,nn) = Binning( proj(:,:,nn), phase_bin ) / phase_bin^2;
    end
    proj = proj_bin;
    clear proj_bin;
    tomo.rot_axis.position = tomo.rot_axis.position / phase_bin;
    tomo.rot_axis.offset = tomo.rot_axis.offset / phase_bin;
    [tomo.vol_shape, tomo.vol_size] = volshape_volsize( proj, [], [], 0);
    prnt( ' done in %g s (%.2f min)', toc - t, (toc - t) / 60)
end

%% Tomographic reco %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tomo.run
    prnt( '\nTomographic reconstruction:')
    
    vol = zeros( tomo.vol_shape, 'single' );

    prnt( '\n volume shape: [%g, %g, %g]', tomo.vol_shape )
    prnt( '\n volume memory allocated: %.2f GiB', Bytes( vol, 3 ) )
    %prod( tomo.vol_shape ) * 4 / 1024^3 ;
    
    if stitch_projections(1)
        rot_axis_offset_reco = 0;
    elseif crop_at_rot_axis(1)
        rot_axis_offset_reco = tomo.rot_axis.position - size( proj, 1) / 2;
    else
        rot_axis_offset_reco = tomo.rot_axis.offset;
    end
    
    if isempty( tomo.take_neg_log )
        tomo.take_neg_log = 1;
    end
    
    % Delete redundant projection
    angles_reco = angles;
    if isequal( angles(1), angles(end) )
        angles_reco(end) = [];
        proj(:,:,end) = [];
    end
    
    % Filter sinogram
    if strcmpi( tomo.algorithm, 'fbp' )
        %%% Provide butterworth filtering also for iterative reconstruction !!!!
        prnt( '\n Filter sino:' )
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
        prnt( ' done in %.2f min.', (toc - t2) / 60)
    end
    
    if crop_at_rot_axis(1)
        % half weight pixel at rot axis pos as it is used twice
        switch excentric_rot_axis
            case 1
                proj( end, :, :) = 0.5 * proj( end, :, :) ;
            case -1
                proj( 1, :, :) = 0.5 * proj( 1, :, :) ;
        end
    end
    
    % Backprojection
    switch lower( tomo.reco_mode )
        case '3d'
            prnt( '\n Backproject:')
            t2 = toc;
            %%% Move parameters to appropiate position or replace globally !!!!!!!!!!!!!!!!!!
            tomo.tilt_camera = ~interactive_mode.lamino * tomo.rot_axis.tilt;
            tomo.tilt_lamino = interactive_mode.lamino * tomo.rot_axis.tilt;
            tomo.angles = tomo.rot_angle.offset + angles_reco;
            %%%%%%%%% Change offset_shit additon %%%%%%%%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tomo.rot_axis.offset = rot_axis_offset_reco + offset_shift;
            vol = astra_parallel3D( tomo, permute(proj, [1 3 2]) );
            pause(0.01)
            prnt( ' done in %.2f min.', (toc - t2) / 60)
            
            vol_min = min( vol(:) );
            vol_max = max( vol(:) );
            
            % Show orthogonal vol cuts
            if visual_output(1)
                
                figure('Name', 'Volume cut z');
                nn = round( size( vol, 3 ) / 2);
                imsc( squeeze( vol(:,:,nn) ) )
                axis equal tight
                title( sprintf( 'vol z = %u', nn ) )
                colorbar
                
                figure('Name', 'Volume cut y');
                nn = round( size( vol, 2 ) / 2);
                im = squeeze( vol(:,nn,:) );
                if size( im , 1) < size( im , 2)
                    imsc( im )
                else
                    imsc( im' )
                end
                axis equal tight
                title( sprintf( 'vol y = %u', nn ) )
                colorbar
                
                figure('Name', 'Volume cut x');
                nn = round( size( vol, 1 ) / 2);
                im = squeeze( vol(nn,:,:) );
                if size( im , 1) < size( im , 2)
                    imsc( im )
                else
                    imsc( im' )
                end
                axis equal tight
                title( sprintf( 'vol x = %u', nn ) )
                colorbar
                
                drawnow
            end
            
            if phase_retrieval.apply
                reco_path = reco_phase_path;
            end
            
            % Save ortho slices x
            nn = round( size( vol, 1 ) / 2);
            im = squeeze( vol(nn,:,:) );
            filename = sprintf( '%s/reco_x%05u', reco_path, nn);
            write32bitTIFfromSingle( filename, im);
            
            % Save ortho slices y
            nn = round( size( vol, 2 ) / 2);
            im = squeeze( vol(:,nn,:) );
            filename = sprintf( '%s/reco_y%05u', reco_path, nn);
            write32bitTIFfromSingle( filename, im);
            
            % Save volume
            if write.reco
                reco_bin = write.reco_binning_factor; % alias for readablity
                CheckAndMakePath( reco_path )
                
                % Save reco path to file
                filename = [userpath, filesep, 'experiments/p05/path_to_latest_reco'];
                fid = fopen( filename , 'w' );
                fprintf( fid, '%s', reco_path );
                fclose( fid );
                
                % Single precision: 32-bit float tiff
                write_volume( write.float, vol, 'float', reco_path, raw_bin, phase_bin, 1, 0, verbose);
                
                % Compression of dynamic range
                if write.uint8 || write.uint8_binned || write.uint16 || write.uint16_binned
                    [tlow, thigh] = compression( vol, write.compression.method, write.compression.parameter );
                else
                    tlow = 0;
                    thigh = 1;
                end
                
                % 16-bit tiff
                write_volume( write.uint16, (vol - tlow)/(thigh - tlow), 'uint16', reco_path, raw_bin, phase_bin, 1, 0, verbose);
                
                % 8-bit tiff
                write_volume( write.uint8, (vol - tlow)/(thigh - tlow), 'uint8', reco_path, raw_bin, phase_bin, 1, 0, verbose);
                
                % Bin data
                if write.float_binned || write.uint16_binned || write.uint8_binned || write.uint8_segmented
                    prnt( '\n Binning:')
                    t2 = toc;
                    vol = Binning( vol, reco_bin ) / reco_bin^3;
                    prnt( ' done in %.2f min.', (toc - t2) / 60)
                end
                
                % Binned single precision: 32-bit float tiff
                write_volume( write.float_binned, vol, 'float', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose);
                
                % 16-bit tiff binned
                write_volume( write.uint16_binned, (vol - tlow)/(thigh - tlow), 'uint16', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose);
                
                % 8-bit tiff binned
                write_volume( write.uint8_binned, (vol - tlow)/(thigh - tlow), 'uint8', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose);
                
                % segmentation
                if write.uint8_segmented
                    [vol, out] = segment_volume(vol, 2^10, visual_output, verbose);
                    save_path = write_volume( 1, vol/255, 'uint8', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose, '_segmented');
                    save( sprintf( '%ssegmentation_info.m', save_path), 'out', '-mat', '-v7.3')
                end
                
            end
            prnt( ' done in %.1f s (%.2f min)', toc-t, (toc-t)/60 )
            
        case 'slice'
            %% Slicewise backprojection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            prnt( '\n Backproject and save slices:')
            t2 = toc;
            
            if phase_retrieval.apply
                reco_path = reco_phase_path;
            end
            
            if write.reco
                reco_bin = write.reco_binning_factor; % alias for readablity
                CheckAndMakePath( reco_path )
                % Save reco path to file
                filename = [userpath, filesep, 'experiments/p05/path_to_latest_reco'];
                fid = fopen( filename , 'w' );
                fprintf( fid, '%s', reco_path );
                fclose( fid );
            end
            
            tomo.angles = tomo.rot_angle.offset + angles_reco;
            %%%%%%%%% Change offset_shit additon %%%%%%%%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tomo.rot_axis.offset = rot_axis_offset_reco + offset_shift;
            
            % Reconstruct central slices first
            [~, ind] = sort( abs( (1:size( proj, 2)) - round( size( proj, 2) / 2 ) ) );
            
            % Loop over slices
            vol_min = Inf;
            vol_max = -Inf;
            for nn = 1:size( proj, 2 )
                
                % Backproject
                vol = astra_parallel2D( tomo, permute( proj(:,ind(nn),:), [3 1 2]) );
                
                vol_min = min( vol_min, min( vol(:) ) );
                vol_max = max( vol_max, max( vol(:) ) );
                
                % Show orthogonal vol cuts
                if visual_output(1) && ~mod( nn - 1, 100 )
                    figure('Name', sprintf( 'Volume slice'));
                    imsc( vol )
                    axis equal tight
                    title( sprintf( 'vol z = %u', ind(1) ) )
                    colorbar
                    drawnow
                    pause( 0.1 )
                end
                
                % Save
                %%% ADD other format options!!!!!!!!!!!!!!!!!!!!!!!!!
                if write.reco
                    % Single precision: 32-bit float tiff
                    write_volume( write.float, vol, 'float', reco_path, raw_bin, phase_bin, 1, ind(nn) - 1, 0);
                end
                
            end
            prnt( ' done in %.2f min.', (toc - t2) / 60)
            
    end
end

%% Write reco log file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if write.reco
    logfile_path = reco_path;
else
    logfile_path = write.path;
end
if write.reco
    logfile_name = sprintf( '%sreco.log', logfile_path);
    fid = fopen( logfile_name, 'w');
    fprintf(fid, 'scan_name : %s\n', scan_name);
    fprintf(fid, 'beamtime_id : %s\n', beamtime_id);
    fprintf(fid, 'scan_path : %s\n', scan_path);
    fprintf(fid, 'reco_path : %s\n', reco_path);
    fprintf(fid, 'MATLAB notation, index of first element: 1, range: first:stride:last\n');
    fprintf(fid, 'MATLAB version : %s\n', version);
    fprintf(fid, 'platform : %s\n', computer);
    fprintf(fid, 'camera : %s\n', cam);
    fprintf(fid, 'num_dark_found : %u\n', num_dark);
    fprintf(fid, 'num_ref_found : %u\n', num_ref_found);
    fprintf(fid, 'num_ref_used : %u\n', num_ref_used);
    fprintf(fid, 'ref_range : %u:%u:%u\n', ref_range(1), ref_range(2) - ref_range(1), ref_range(end) );
    fprintf(fid, 'num_proj_found : %u\n', num_proj_found);
    fprintf(fid, 'num_proj_used : %u\n', num_proj_used);
    fprintf(fid, 'proj_range : %u:%u:%u\n', proj_range(1), proj_range(2) - proj_range(1), proj_range(end) );
    fprintf(fid, 'im_shape_raw : %u %u\n', im_shape_raw);
    fprintf(fid, 'raw_roi : %u %u\n', raw_roi);
    fprintf(fid, 'image_shape_roi : %u %u\n', size( im_roi ));
    fprintf(fid, 'image_shape_roi_binned : %u %u\n', im_shape_binned);
    fprintf(fid, 'raw_binning_factor : %u\n', raw_bin);
    fprintf(fid, 'bin_before_filtering : %u\n', bin_before_filtering);
    fprintf(fid, 'effective_pixel_size : %g micron\n', eff_pixel_size * 1e6);
    fprintf(fid, 'effective_pixel_size_binned : %g micron\n', eff_pixel_size_binned * 1e6);
    fprintf(fid, 'energy : %g eV\n', energy);
    fprintf(fid, 'sample_detector_distance : %f m\n', sample_detector_distance);
    fprintf(fid, 'ring_current_normalization : %u\n', ring_current_normalization);
    fprintf(fid, 'image_correlation.method : %s\n', image_correlation.method);
    fprintf(fid, 'image_correlation.num_flats : %u\n', image_correlation.num_flats);
    fprintf(fid, 'image_correlation.area_width : %u:%u:%u\n', image_correlation.area_width(1), image_correlation.area_width(2) - image_correlation.area_width(1), image_correlation.area_width(end));
    fprintf(fid, 'image_correlation.area_height : %u:%u:%u\n', image_correlation.area_height(1), image_correlation.area_height(2) - image_correlation.area_height(1), image_correlation.area_height(end));
    if ~read_flatcor
        fprintf(fid, 'min_max_of_all_darks : %6g %6g\n', darks_min, darks_max);
        fprintf(fid, 'min_max_of_median_dark : %6g %6g\n', dark_med_min, dark_med_max);
        fprintf(fid, 'min_max_of_all_flats : %6g %6g\n', flat_min, flat_max);
        fprintf(fid, 'min_max_of_all_corrected_flats : %6g %6g\n', flat_min2, flat_max2);
        fprintf(fid, 'min_max_of_all_raws :  %6g %6g\n', raw_min, raw_max);
        fprintf(fid, 'min_max_of_all_corrected_raws :  %6g %6g\n', raw_min2, raw_max2);
        fprintf(fid, 'min_max_of_all_flat_corr_projs : %g %g \n', proj_min, proj_max);
    end
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
    % Volume
    fprintf(fid, 'tomo.vol_shape : %u %u %u\n', tomo.vol_shape(1), tomo.vol_shape(2), tomo.vol_shape(3));
    fprintf(fid, 'tomo.vol_size : %f %f %f %f %f %f\n', tomo.vol_size(1), tomo.vol_size(2), tomo.vol_size(3), tomo.vol_size(4), tomo.vol_size(5), tomo.vol_size(6));
    % Rotation
    fprintf(fid, 'excentric_rot_axis : %i\n', excentric_rot_axis);
    fprintf(fid, 'crop_at_rot_axis : %u\n', crop_at_rot_axis);
    fprintf(fid, 'stitch_projections : %u\n', stitch_projections);
    fprintf(fid, 'stitch_method : %s\n', stitch_method );
    fprintf(fid, 'tomo.reco_mode : %s\n', tomo.reco_mode);
    fprintf(fid, 'tomo.rot_angle.full_range : %f * pi rad\n', tomo.rot_angle.full_range / pi);
    fprintf(fid, 'tomo.rot_angle.offset : %f * pi rad\n', tomo.rot_angle.offset / pi);
    fprintf(fid, 'rotation_axis_offset_calculated : %f\n', rot_axis_offset_calc);
    fprintf(fid, 'rot_axis_offset_reco : %f\n', rot_axis_offset_reco);
    fprintf(fid, 'rotation_axis_offset_reco : %f\n', rot_axis_offset_reco);
    fprintf(fid, 'rotation_axis_position_calculated : %f\n', rot_axis_pos_calc);
    fprintf(fid, 'tomo.rot_axis.position : %f\n', tomo.rot_axis.position);
    fprintf(fid, 'rotation_axis_tilt_calculated : %f\n', rot_axis_tilt_calc);
    fprintf(fid, 'tomo.rot_axis.tilt : %f\n', tomo.rot_axis.tilt);
    fprintf(fid, 'raw_image_binned_center : %f\n', im_shape_binned1 / 2);
    fprintf(fid, 'tomo.rot_axis.corr_area1 : %u:%u:%u\n', tomo.rot_axis.corr_area1(1), tomo.rot_axis.corr_area1(2) - tomo.rot_axis.corr_area1(1), tomo.rot_axis.corr_area1(end));
    fprintf(fid, 'tomo.rot_axis.corr_area2 : %u:%u:%u\n', tomo.rot_axis.corr_area2(1), tomo.rot_axis.corr_area2(2) - tomo.rot_axis.corr_area2(1), tomo.rot_axis.corr_area2(end));
    fprintf(fid, 'interactive_mode.rot_axis_pos : %u\n', interactive_mode.rot_axis_pos);
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
    % FBP
    fprintf(fid, 'tomo.fbp_filter.type : %s\n', tomo.fbp_filter.type);
    fprintf(fid, 'tomo.fbp_filter.freq_cutoff : %f\n', tomo.fbp_filter.freq_cutoff);
    fprintf(fid, 'tomo.fbp_filter.padding : %u\n', tomo.fbp_filter.padding);
    fprintf(fid, 'tomo.fbp_filter.padding_method : %s\n', tomo.fbp_filter.padding_method);
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
        fprintf(fid, 'write.compression.method : %s\n', write.compression.method);
        if exist( 'tlow', 'var' ) && exist( 'thigh', 'var' )
            fprintf(fid, 'compression_limits : %f %f\n', tlow, thigh);
        end
        fprintf(fid, 'reco_bin : %u\n', reco_bin);
    end
    fprintf(fid, 'full_reconstruction_time : %.1f s\n', toc);
    fprintf(fid, 'date_of_reconstruction : %s\n', datetime);
    fprintf(fid, 'tomo.rot_axis.offset at %u x binning : %f\n', raw_bin, rot_axis_offset_reco);
    fclose(fid);
    % End of log file
    prnt( '\n log file : %s', logfile_name)
    prnt( '\n reco_path : \n%s', reco_path)
end
PrintVerbose(verbose && interactive_mode.rot_axis_pos, '\nTime elapsed in interactive mode: %g s (%.2f min)', tint, tint / 60 );
prnt( '\nTime elapsed for computation: %g s (%.2f min)', toc - tint, (toc - tint) / 60 );
prnt( '\nFINISHED: %s\n\n', scan_name)
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dbclear if error
