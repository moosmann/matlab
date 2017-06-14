% P05 reconstruction pipeline: preprocessing, filtering, phase retrieval,
% tomographic reconstruction, etc.
%
% USAGE
% Set parameters in PARAMETERS / SETTINGS section below and run script.
%
% How to run script:
% - type 'F5' when focus is in the Editor windows
% - click 'Run' in the Editor tab
% - type 'p05_reco_loop' and Enter in the Command Window
%
% To loop over sets of data or parameters sets use 'p05_reco_loop'.
%
% Written by Julian Moosmann. First version: 2016-09-28. Last modifcation:
% 2017-06-06

close all hidden % open windows
%dbstop if error

%% PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scan_path = ...    
    '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/syn13_55L_Mg10Gd_12w_load_00';        
    '/asap3/petra3/gpfs/p05/2017/data/11002839/raw/ehh_2017_015_a';
'/asap3/petra3/gpfs/p05/2017/data/11003950/raw/syn01_48L_PEEK_12w_b';    
read_flatcor = 0; % read flatfield-corrected images from disc, skips preprocessing
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = [];[201 2400]; % [y0 y1] vertical roi.  skips first raw_roi(1)-1 lines, reads until raw_roi(2)
raw_bin = 2; % projection binning factor: 1, 2, or 4
excentric_rot_axis = 0; % off-centered rotation axis increasing FOV. -1: left, 0: centeerd, 1: right. influences rot_corr_area1
crop_at_rot_axis = 0; % recommended for scans with excentric rotation axis when no projection stitching is done
stitch_projections = 0; % stitch projection (for 2 pi scans) at rotation axis position. "doubles" number of voxels
stitch_method = 'sine'; % 'step': no interpolation, 'linear','sine': linear interpolation of overlap region. !!! adjust: correlation area
proj_range = 1; % range of projections to be used (from all that are found). if empty: all, if scalar: stride
ref_range = 1; % range of flat fields to be used (from all that are found). start:incr:end. if empty: all (equals 1). if scalar: stride
darkFiltPixHot = 0.01; % Hot pixel filter parameter for dark fields, for details see 'FilterPixel'
darkFiltPixDark = 0.005; % Dark pixel filter parameter for dark fields, for details see 'FilterPixel'
refFiltPixHot = 0.01; % Hot pixel filter parameter for flat fields, for details see 'FilterPixel'
refFiltPixDark = 0.005; % Dark pixel filter parameter for flat fields, for details see 'FilterPixel'
projFiltPixHot = 0.01; % Hot pixel filter parameter for projections, for details see 'FilterPixel'
projFiltPixDark = 0.005; % Dark pixel filter parameter for projections, for details see 'FilterPixel'
correlate_proj_flat = 1;%  correlate flat fields and projection to correct beam shaking
correlation_method = 'ssim-ml'; 'diff';'shift';'ssim';'std';'entropy';'cov';'cross'; % method for correlation. 'diff': difference measure (preferred), 'shift': computes relative shift using image cross-correlation
corr_shift_max_pixelshift = 0.25; % maximum pixelshift allowed for 'shift'-correlation method: if 0 use the best match (i.e. the one with the least shift), if > 0 uses all flats with shifts smaller than corr_shift_max_pixelshift
corr_num_flats = 3; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
norm_by_ring_current = 1; % normalize flat fields and projections by ring current
flat_corr_area1 = [1 floor(100/raw_bin)]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
flat_corr_area2 = [0.2 0.8]; %correlation area: index vector or relative/absolute position of [first pix, last pix]
decimal_round_precision = 2; % precision when rounding pixel shifts
ring_filter = 1; % ring artifact filter
ring_filter_method = 'jm';'wavelet-fft'; 
ring_filter_median_width = 11; %[3 11 21 31 39];
dec_levels = 2:6;
wname = 'db30';'db25';
sigma = 2.4;
% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_phase_retrieval = 0;
phase_retrieval_method = 'tie';'qp';'qpcut'; 'tie'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval_reg_par = 2.5; % regularization parameter
phase_retrieval_bin_filt = 0.15; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval_cutoff_frequ = 1 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_padding = 1; % padding of intensities before phase retrieval
energy = []; % in eV. if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
% Tomography parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_tomo = 1; % run tomographic reconstruction
vol_shape = [];% shape of the volume to be reconstructed, either in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
vol_size = []; % if empty, unit voxel size is assumed
rot_angle_full = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
rot_angle_offset = pi; % global rotation of reconstructed volume
rot_axis_offset = [] ; % if empty use automatic computation
rot_axis_pos = []; % if empty use automatic computation. either offset or pos has to be empty. can't use both
rot_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axisq
rot_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
rot_corr_gradient = 0; % use gradient of intensity maps if signal variations are too weak to correlate projections
rot_axis_tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
fbp_filter_type = 'Ram-Lak';'linear';
fpb_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions
fbp_filter_padding_method = 'symmetric';
butterworth_filter = 1; % use butterworth filter in addition to FBP filter
butterworth_order = 1;
butterworth_cutoff_frequ = 0.5;
astra_pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted
take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '';% absolute path were output data will be stored. !!overwrites the write_to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write_to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write_flatcor = 0; % save preprocessed flat corrected projections
write_phase_map = 0; % save phase maps (if phase retrieval is not 0)
write_sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write_sino_phase = 0; % save sinograms of phase mapsls
write_reco = 1; % save reconstructed slices (if do_tomo=1)
write_float = 0; % write single precision (32-bit float) tiff
write_float_binned = 1; % write binned single precision (32-bit float) tiff
write_16bit = 0; % write 8bit-tiff, currently for reco only
write_8bit = 0; % write 8bit-tiff, currently for reco only
write_8bit_binned = 0; % write binned 8bit-tiff, currently for reco only
reco_bin = 2; % currently only 2 x 2 x 2 binning is implemented
compression = 'histo';'full'; 'std'; 'threshold'; % method of compression of dynamic range
compression_std_num = 5; % dynamic range: mean(volume) +/- NUM* std(volume)
compression_threshold = [-1 1]; % dynamic range: [MIN MAX]
compression_histo = [0.02 0.02]; % [LOW HIGH]. crop dynamic range to values between (100*LOW)% and (100*HIGH)% of the original histogram
parfolder = ''; % parent folder for 'reco', 'sino', 'phase', and 'flat_corrected'
subfolder_flatcor = ''; % subfolder in 'flat_corrected'
subfolder_phase_map = ''; % subfolder in 'phase_map'
subfolder_sino = ''; % subfolder in 'sino'
subfolder_reco = ''; % subfolder in 'reco'
% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visualOutput = 1; % show images and plots during reconstruction
interactive_determination_of_rot_axis = 1; % reconstruct slices with different rotation axis offsets
interactive_determination_of_rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_determination_of_rot_axis_slice = 0.5; % slice number, default: 0.5. if in [0,1): relative, if in (1, N]: absolute
% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
poolsize = 0.80; % number of workers used in a parallel pool. if > 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used
link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all

%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: physically consistent attenutation values of reconstructed slice 
% TODO: check raw roi, transpose, rot90 for KIT and EHD camera
% TODO: get rid of transpose when reading image files
% TODO: vertical stitching
% TODO: volume shape for excentric rot axis
% TODO: reco bin once
% TODO: Interactive plot mode for wavelet-fft fing filter
% TODO: Improve FilterStripesCombinedWaveletFFT
% TODO: check proj-flat correlation measure for absolute values
% TODO: check crop_at_rot_axis option with stitch_projections, etc
% TODO: SSIM: include Gaussian blur filter, test
% TODO: large data set management: parloop, memory, etc
% TODO: stitching: optimize and refactor, memory efficency, interpolation method
% TODO: interactive loop over tomo slices for different phase retrieval parameter
% TODO: automatic determination of rot center
% TODO: output file format option: 8-bit, 16-bit for all saved images
% TODO: additional padding schemes for FBP filter
% TODO: read sinogram option
% TODO: set photometric tag for tif files w/o one, turn on respective warning
% TODO: GPU phase retrieval: parfor-loop requires memory managment
% TODO: median filter width of ring filter dependence on binning
% TODO: check offset: proj correlation for rotation axis determination
% TODO: check offset: flat/proj correlation 
% TODO: inverse Gaussian filter for phase retrieval, VZC theorem
% TODO: flat-flat correlation and averaging before proj-flat correlation
% TODO: proj-flat correlation: take negative logarithm or not?

%% Notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Padding and phase retrieval:
% Symmetric padding of intensity maps before phase retrieval cleary reduces
% artifacts in the retrieved phase maps which are due to inconistent
% boundary conditions. However, due to local tomography phase retrieval
% with padding can result in more artifacts in the recostructed volume
% (such as additonal bumps at the halo described below). After phase
% retrieval without padding the region close to the left
% boundary is blended in to the region close to the right boundary and vice
% versa (same holds for top and bottom). This effect fades out with
% increasing distance to the boundarys. Due to this blending effect
% asymmetrical 'density' distributions, which are the cause of halo-like
% artifacts stretching outwards from the biggest possible circle within the
% reconstruction volume within a slice, are reduced which reduces the
% halo-like artifacts as well.

% FOV extension by excentric rotation axis:
% For absorpion-contrast data or more precisely when no phase retrieval
% is desired, volumes can be reconstructed from a data set where an excentric
% position of the rotation axis was used without prior stitching of the
% projections by simply providing the correct rotation axis position and
% setting fbp_filter_padding to 1. Maybe, for the automatic detection of the
% rotation axis to work, the area to correlate has to be adjusted i.e.
% 'rot_corr_area1'. When using phase retrieval, this
% approach does not work appropriately and gives rise to artifacts near the
% center of the reconstructed volume. This is due to the fact, that without
% stitching the phase is retrieved from a 'cropped' projection which
% results in inconsitently retrieved low frequencies (large scale
% variations) in the phase map. Using the 'linear' FBP filter instead of
% 'Ram-Lak' maybe reduces these artifacts (not tested).

% Correlation of projections and flat fields:
% Compared to the simple differencing method, the cross correlation method
% is very likely to give non-optimal results since is more sensitive to
% small-scale features, such as those stemming from contimation of the beam
% from the scintllator, diamond window, etc, and less sensitive to
% variations on a larger scale.

% Entropy-type determination of rotation axis position
% Empirically, this does not work well for excentric rotation axis positions and
% stitchted projections. To be tested

% Ring artefact filter

% GPU CUDA error: device busy
% Do not call 'gpuDevice()' when ASTRA toolbox is used with multi-GPU
% support, or MATLAB terminates abnormally with a segmentation violation.

%% External call: parameters set by 'p05_reco_loop' %%%%%%%%%%%%%%%%%%%%%%%
if exist( 'external_parameter' ,'var')
    visualOutput = 0;
    interactive_determination_of_rot_axis = 0;
    fields = fieldnames( external_parameter );
    for nn = 1:numel( fields )
        var_name = fields{nn};
        var_val = getfield( external_parameter, var_name );
        assignin('caller', var_name, var_val )
    end
    clear external_parameter;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
PrintVerbose(verbose, 'Start reconstruction of ')
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
cdscandir = cd( scan_path );
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';
imsc1 = @(im) imsc( flipud( im' ) );
astra_clear % if reco was aborted, ASTRA memory is not cleared
if ~isempty( rot_axis_offset ) && ~isempty( rot_axis_pos )
    error('rot_axis_offset (%f) and rot_axis_pos (%f) cannot be used simultaneously. One must be empty.', rot_axis_offset, rot_axis_pos)
end

%% Input folder
while scan_path(end) == filesep
    scan_path(end) = [];
end
[raw_path, scan_name] = fileparts(scan_path);
% Save raw path in file
filename = [userpath, filesep, 'experiments/p05/pathtolastraw'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', raw_path );
fclose( fid );
% Scan path
scan_path = [scan_path, filesep];
[beamtime_path, raw_folder] = fileparts(raw_path);
[~, beamtime_id] = fileparts(beamtime_path);
if ~strcmp(raw_folder, 'raw')
    error('Name of folder is not raw: %s', raw_folder)
end
PrintVerbose(verbose, '%s', scan_name)
PrintVerbose(verbose, '\n scan_path:%s', scan_path)
% Save scan path in file
filename = [userpath, filesep, 'experiments/p05/pathtolastscan'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', scan_path );
fclose( fid );

%% Output folder
out_folder = 'processed';
if write_to_scratch(1)
    out_folder = 'scratch_cc';
end
if isempty( out_path )
    out_path = [beamtime_path, filesep, out_folder, filesep, scan_name];
else
    out_path = [out_path, filesep, scan_name];
end
if ~isempty(parfolder)
    out_path = [out_path, filesep, parfolder];
end

% Projections path
if isempty( subfolder_flatcor )
    flatcor_path = [out_path, filesep, 'flat_corrected', filesep];
else
    flatcor_path = [out_path, filesep, 'flat_corrected', filesep, subfolder_flatcor, filesep];
end
PrintVerbose(verbose & write_flatcor, '\n flatcor_path: %s', flatcor_path)

% Phase map path
if isempty( subfolder_phase_map )
    phase_map_path = [out_path, filesep, 'phase_map', filesep];
else
    phase_map_path = [out_path, filesep, 'phase_map', filesep, subfolder_phase_map, filesep];
end
PrintVerbose(verbose & write_phase_map, '\n phase_map_path: %s', phase_map_path)

% Sinogram path
if isempty( subfolder_sino )
    sino_path = [out_path, filesep, 'sino', filesep];
    sino_phase_path = [out_path, filesep, 'sino_phase', filesep];
else
    sino_path = [out_path, filesep, 'sino', filesep, subfolder_sino, filesep];
    sino_phase_path = [out_path, filesep, 'sino_phase', filesep, subfolder_sino, filesep];
end
PrintVerbose(verbose & write_sino, '\n sino_path: %s', sino_path)
PrintVerbose(verbose & write_sino_phase, '\n sino_phase_path: %s', sino_phase_path)

% Reco path
if isempty( subfolder_reco )
    reco_path = [out_path, filesep, 'reco', filesep];
else
    reco_path = [out_path, filesep, 'reco', filesep, subfolder_reco, filesep];
end
PrintVerbose(verbose, '\n reco_path : %s', reco_path)

% Memory
[mem_free, mem_avail, mem_total] = free_memory;
PrintVerbose(verbose, '\n system memory: free, available, total : %.3g GiB, %.3g GiB, %.3g GiB', mem_free/1024^3, mem_avail/1024^3, mem_total/1024^3)
if isempty( gpu_index )
    gpu_index = 1:gpuDeviceCount;
end
for nn = gpu_index
    gpu = parallel.gpu.GPUDevice.getDevice( nn );
    PrintVerbose(verbose, '\n gpu %u : memory: available, total, percent : %.3g GiB, %.3g GiB, %.2f%%', nn, gpu.AvailableMemory/1024^3, gpu.TotalMemory/1024^3, 100*gpu.AvailableMemory/gpu.TotalMemory)
end

%% File names

% Projection file names
proj_names = FilenameCell( [scan_path, '*.img'] );
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, 'proj_*.tif'] );
end
num_proj_found = numel(proj_names);

% Ref file names
ref_names = FilenameCell( [scan_path, '*.ref'] );
if isempty( ref_names )
    ref_names = FilenameCell( [scan_path, 'ref_*.tif'] );
end
num_ref_found = numel(ref_names);
if isempty( ref_range )
    ref_range = 1;
end
if numel( ref_range ) == 1
    ref_range = 1:ref_range:num_ref_found;
end
num_ref_used = numel( ref_range );
ref_names_mat = NameCellToMat( ref_names(ref_range) );
ref_nums = CellString2Vec( ref_names(ref_range) );
PrintVerbose(verbose, '\n number of refs found : %g', num_ref_found)
PrintVerbose(verbose, '\n number of refs used : %g', num_ref_used)
PrintVerbose(verbose, '\n reference range used : %g:%g:%g%', ref_range(1), ref_range(2) - ref_range(1), ref_range(end))

% Dark file names
dark_names = FilenameCell( [scan_path, '*.dar'] );
if isempty( dark_names )
    dark_names = FilenameCell( [scan_path, 'dark_*.tif'] );
end
dark_nums = CellString2Vec( dark_names );
num_dark = numel(dark_names);
PrintVerbose(verbose, '\n number of darks found : %g', num_dark)

% Projection range to read
if isempty( proj_range )
    proj_range = 1;
end
if numel( proj_range ) == 1
    proj_range = 1:proj_range:num_proj_found;
end
num_proj_used = numel( proj_range );
proj_nums = CellString2Vec( proj_names(proj_range) );
PrintVerbose(verbose, '\n number of projections found : %g', num_proj_found)
PrintVerbose(verbose, '\n number of projections used : %g', num_proj_used)
PrintVerbose(verbose, '\n projection range used : first:stride:last =  %g:%g:%g', proj_range(1), proj_range(2) - proj_range(1), proj_range(end))

%% Image shape and ROI
filename = sprintf('%s%s', scan_path, ref_names{1});
[im_raw, tif_info] = read_image( filename );
raw_im_shape_raw = size( im_raw );
im_roi = read_image( filename, '', raw_roi, tif_info );
raw_im_shape = size( im_roi );
raw_im_shape_binned = floor( raw_im_shape / raw_bin );
raw_im_shape_binned1 = raw_im_shape_binned(1);
raw_im_shape_binned2 = raw_im_shape_binned(2);
PrintVerbose(verbose, '\n raw image shape : %g  %g', raw_im_shape_raw)
PrintVerbose(verbose, '\n raw image shape roi : %g  %g', raw_im_shape)
PrintVerbose(verbose, '\n raw image shape binned : %g  %g', raw_im_shape_binned)

%% P05 log-file
str = dir( sprintf( '%s*scan.log', scan_path) );
filename = sprintf( '%s/%s', str.folder, str.name);
[par, cur, cam] = p05_log( filename );
if isempty( energy )
    if isfield( par, 'Energy' )
        energy = par.Energy;
    elseif isfield( par, 'energy' )
        energy = par.energy;
    end
end
if isempty( eff_pixel_size )
    if isfield( par, 'eff_pix_size' )
        eff_pixel_size = par.eff_pix_size * 1e-3 ;
    elseif isfield( par, 'eff_pix' )
        eff_pixel_size = par.eff_pix * 1e-3 ;
    elseif isfield( par, 'ccd_pixsize' ) && isfield( par, 'magn' )
        eff_pixel_size = par.ccd_pixsize / par.magn * 1e-3 ;
    end
end
eff_pixel_size_binned = raw_bin * eff_pixel_size;
if isempty( sample_detector_distance )
    if isfield( par, 'camera_distance')
        sample_detector_distance = par.camera_distance / 1000;
    elseif isfield( par, 'o_ccd_dist')
        sample_detector_distance = par.o_ccd_dist / 1000;
    end
end
PrintVerbose( verbose, '\n energy : %.1f keV', energy / 1e3 )
PrintVerbose( verbose, '\n distance sample dector : %.1f mm', sample_detector_distance * 1000 )
PrintVerbose( verbose, '\n effective pixel size : %.2f micron',  eff_pixel_size * 1e6)

%% Start parallel CPU pool
t = toc;
if poolsize <= 1 && poolsize > 0
    poolsize = max( floor( poolsize * feature('numCores') ), 1 );
end
PrintVerbose( verbose && (poolsize > 1), '\nStart parallel pool of %u workers. ', poolsize)
OpenParpool(poolsize);
PrintVerbose( verbose && (poolsize > 1), ' Time elapsed: %.1f s', toc-t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read flat corrected projection
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
    PrintVerbose(verbose, '\n Read flat corrected projections.')
    if num_proj_read ~= num_proj_used
        fprintf('\n WARNING: Number of flat corrected projections read (%g) differs from number of projections to be processed (%g)!\n', num_proj_read, num_proj_used)
    end
    
    % Preallocation
    proj = zeros( raw_im_shape_binned(1), raw_im_shape_binned(2), num_proj_read, 'single');
    PrintVerbose( verbose, ' Allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
    proj_names_mat = NameCellToMat( proj_names );
    
    % Read flat corrected projections
    parfor nn = 1:num_proj_read
        filename = sprintf('%s%s', flatcor_path, proj_names_mat(nn, :));
        proj(:, :, nn) = fliplr( read_image( filename ) );
    end
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    
    %% Read raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~read_flatcor
    %% Dark field
    t = toc;
    PrintVerbose(verbose, '\nProcessing dark fields.')
    darks = zeros( [raw_im_shape_binned, num_dark], 'single');
    PrintVerbose( verbose, ' Allocated bytes: %.2f MiB.', Bytes( darks, 2 ) )
    parfor nn = 1:num_dark
        filename = sprintf('%s%s', scan_path, dark_names{nn});
        im = single( read_image( filename, '', raw_roi, tif_info) );
        % Remove large outliers. Assume Poisson distribtion at large lambda
        % is approximately a Gaussian distribution and set all value above
        % mean + 4 * std (99.994 of values lie within 4 std). Due to
        % outliers 4*std will contain much more values and is a good
        % estimate
        im_mean = mean( im(:) );
        im_std = std( im(:) );
        im( im > im_mean + 4*im_std) = im_mean;
        darks(:, :, nn) = Binning( FilterPixel( im, [darkFiltPixHot darkFiltPixDark]), raw_bin) / raw_bin^2;
    end
    dark_min = min( darks(:) );
    dark_max = max( darks(:) );
    dark = squeeze( median(darks, 3) );
    dark_med_min = min( dark(:) );
    dark_med_max = max( dark(:) );
    PrintVerbose(verbose, ' Time elapsed: %.1f s', toc-t)
    if visualOutput(1)
        h1 = figure('Name', 'mean dark field, flat field, projections');
        subplot(2,2,1)
        imsc1( dark );
        axis equal tight
        title(sprintf('dark field'))
        colorbar
        drawnow
    end
    
    %% Flat field
    t = toc;
    PrintVerbose(verbose, '\nProcessing flat fields.')
    
    % Preallocation
    flat = zeros( [raw_im_shape_binned, num_ref_used], 'single');
    PrintVerbose( verbose, ' Allocated bytes: %.2f GiB.', Bytes( flat, 3 ) )
    
    % Parallel loop
    parfor nn = 1:num_ref_used
        filename = sprintf('%s%s', scan_path, ref_names_mat(nn, :));
        flat(:, :, nn) = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info ), [refFiltPixHot refFiltPixDark]), raw_bin) / raw_bin^2;
        % Check for zeros
        num_zeros =  sum( sum( flat(:,:,nn) < 1 ) );
        if num_zeros > 0
            fprintf('\n WARNING: values of %u pixels of flat field no %u are lower than 1.', num_zeros, nn)
        end
    end
    % min/max values before dark field subtraction and ring normalization
    flat_min = min( flat(:) );
    flat_max = max( flat(:) );
    
    % Dark field correction
    flat = bsxfun( @minus, flat, dark );
    
    % Ring current normalization
    if norm_by_ring_current(1)
        switch lower( cam )
            case 'kit'
                ref_ind = ref_nums + 1;
            case 'ehd'
                ref_ind = ref_range;
        end
        if isequal( ref_nums, [cur.ref(ref_ind).ind])
            scale_factor = 100 ./ shiftdim( [cur.ref(ref_ind).val], -1 );
            flat = bsxfun( @times, flat, scale_factor );
            if visualOutput(1)
                hrc = figure('Name', 'Ring current normalization');
                subplot(2,1,1);
                plot( scale_factor(:), '.' )
                axis tight% equal tight square
                title(sprintf('flats'))
                drawnow
            end
        else
            fprintf('\n WARNING: flat fields not normalized by ring current. Names read from dir() and log-file are inconsistent.')
        end
    end
    
    % Check for non-positive values
    parfor nn = 1:num_ref_used
        im = flat(:,:,nn);
        m = im < 1;
        im(m) = 1;
        flat(:,:,nn) = im;
    end
    flat_min2 = min( flat(:) );
    flat_max2 = max( flat(:) );
    
    num_zeros =  sum( flat(:) < 1 );
    if num_zeros > 0
        fprintf('\n WARNING: flat field contains %u zeros', num_zeros)
    end
    
    PrintVerbose(verbose, ' Time elapsed: %.1f s', toc-t)
    if visualOutput(1)
        figure(h1)
        subplot(2,2,2)
        imsc1( flat(:,:,1) )
        axis equal tight% square
        title(sprintf('flat field #1'))
        drawnow
        colorbar
    end
    
    %% Projections
    t = toc;
    PrintVerbose(verbose, '\nRead and filter projections.')
    % Preallocation
    proj = zeros( raw_im_shape_binned(1), raw_im_shape_binned(2), num_proj_used, 'single');
    PrintVerbose( verbose, ' Allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
    img_names_mat = NameCellToMat( proj_names(proj_range) );
    
    % Display first raw images
    if visualOutput(1)
        figure(h1)
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        raw1 = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info ), [projFiltPixHot, projFiltPixDark]), raw_bin) / raw_bin^2;
        subplot(2,2,3)
        imsc1( raw1 )
        axis equal tight
        title(sprintf('raw projection #1 (roi)'))
        drawnow
        colorbar
    end
    
    % Get absolute filter thresholds from percentage-wise pixel filtering
    % of 1st, middle, and last projection to speed up processing
    if projFiltPixHot < 1 || projFiltPixDark < 0.5
        filename = sprintf('%s%s', scan_path, img_names_mat(num_proj_used, :));
        [~, ht(3), dt(3)] = FilterPixel( read_image( filename, '', raw_roi, tif_info ), [projFiltPixHot, projFiltPixDark]);
        
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        [~, ht(2), dt(2)] = FilterPixel( read_image( filename, '', raw_roi, tif_info ), [projFiltPixHot, projFiltPixDark]);
                
        filename = sprintf('%s%s', scan_path, img_names_mat(round(num_proj_used/2), :));
        [~, ht(1), dt(1)] = FilterPixel( read_image( filename, '', raw_roi, tif_info ), [projFiltPixHot, projFiltPixDark]);
        
        HotThresh = mean( ht );
        DarkThresh = mean( dt );        
    else
        HotThresh = projFiltPixHot;
        DarkThresh = projFiltPixDark;                        
    end
    
    % Read raw projections
    parfor nn = 1:num_proj_used
        filename = sprintf('%s%s', scan_path, img_names_mat(nn, :));
        proj(:, :, nn) = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info ), [HotThresh, DarkThresh]), raw_bin) / raw_bin^2;
    end
    raw_min = min( proj(:) );
    raw_max = max( proj(:) );
    
    % Dark field correction
    proj = bsxfun( @minus, proj, dark);
    
    % Ring current normalization
    if norm_by_ring_current(1)
        switch lower( cam )
            case 'kit'
                proj_ind = proj_nums + 1;
            case 'ehd'
                proj_ind = proj_range;
        end
        if isequal( proj_nums, [cur.proj(proj_ind).ind] )
            scale_factor = 100 ./ shiftdim( [cur.proj(proj_ind).val], -1 );
            proj = bsxfun( @times, proj, scale_factor );
            if visualOutput(1)
                figure(hrc)
                subplot(2,1,2);
                plot( scale_factor(:), '.' )
                axis tight %equal
                title(sprintf('projections'))
                drawnow
            end
        else
            fprintf('\n WARNING: projections not normalized by ring current. Names read from dir() and log-file are not consistent.')
        end
    end
    
    % Check for non-positive values
    parfor nn = 1:num_proj_used
        im = proj(:,:,nn);
        m = im < 1;
        if sum( m(:) ) > 0
            im(m) = 1;
            proj(:,:,nn) = im;
        end
    end
    raw_min2 = min( proj(:) );
    raw_max2 = max( proj(:) );
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    PrintVerbose(verbose, '\n hot- / dark-pixel filter threshold : %f, %f', HotThresh, DarkThresh )    
    
    %% Flat field correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if ~correlate_proj_flat
        % Flat field correction without correlation
        flat_m = median(flat, 3);
        proj = bsxfun( @times, proj, flat_m);    
    else

        % Correlate flat fields    
        PrintVerbose(verbose, '\nCorrelate projections and flat-fields.')
        t = toc;
        % Correlation ROI
        flat_corr_area1 = IndexParameterToRange(flat_corr_area1, raw_im_shape_binned(1));
        flat_corr_area2 = IndexParameterToRange(flat_corr_area2, raw_im_shape_binned(2));
        flat_roi = flat(flat_corr_area1, flat_corr_area2, :);
        proj_roi = proj(flat_corr_area1, flat_corr_area2, :);
        flat_corr_shift_1 = zeros( num_proj_used, num_ref_used);
        flat_corr_shift_2 = zeros( num_proj_used, num_ref_used);
        c_diff1_l1 = zeros( num_proj_used, num_ref_used);
        c_diff1_l2 = c_diff1_l1;
        c_diff2_l1 = c_diff1_l1;
        c_diff2_l2 = c_diff1_l1;
        c_std = zeros( num_proj_used, num_ref_used);
        c_ent = c_std;
        c_cov = c_std;
        c_corr = c_std;
        c_ssim = c_std;
        c_ssim_ml = c_std;
        c_l = c_std;
        c_c = c_std;
        c_s = c_std;
        % Compute shift for each pair projection/flat field
        for ff = 1:num_ref_used
            flat_ff = flat_roi(:,:,ff);
            
            switch correlation_method
                case 'shift'
                    parfor pp = 1:num_proj_used
                        % (p^*)(-x) ** f(x): cross correlation equals
                        % convolution of complex conjugate of p(-x), i.e.
                        % rot180(p) and f(x)
                        out = ImageCorrelation(proj_roi(:,:,pp), flat_ff, 0, 0, 0, 0, 1);
                        flat_corr_shift_1(pp,ff) = round( out.shift1, decimal_round_precision );
                        flat_corr_shift_2(pp,ff) = round( out.shift2, decimal_round_precision) ; % relevant shift
                    end
                case {'diff', 'std', 'entropy', 'ssim', 'ssim-ml', 'cov', 'corr'}
                    % Dynamic range of camera for SSIM
%                     switch lower( cam )
%                         case 'kit'
%                             L = 2^12;
%                         case 'ehd'
%                             L = 2^16;
%                     end
                    %L = round(max(max(flat_roi(:)),max(flat_roi(:))) - min(min(flat_roi(:)),min(flat_roi(:))));
                    L = 1;
                    parfor pp = 1:num_proj_used
                        % projection
                        p = proj_roi(:,:,pp); 
                        p_mean = mean2( p ); 
                        p_std = std2( p );
                        % flat field
                        f = flat_ff;
                        f_mean = mean2( f );
                        f_std = std2( f );
                        
                        % differences
                        d1 =  abs( abs( p ) - abs( f ) ) ;
                        d2 =  sqrt( abs( p.^2 - f.^2) );
                        
                        % ratio
                        r =  p ./ f ;
                        
                        % measures
                        % difference: anisotropic L1
                        c_diff1_l1(pp,ff) = norm( d1(:), 1);
                        % difference: anisotropic L2
                        c_diff1_l2(pp,ff) = norm( d1(:), 2);
                        % difference: isotropic L1
                        c_diff2_l1(pp,ff) = norm( d2(:), 1);
                        % difference: isotropic L2
                        c_diff2_l2(pp,ff) = norm( d2(:), 2);                                                
                        % std: standard deviation of ratio
                        c_std(pp,ff) = std2( r );
                        % entropy of ratio
                        c_ent(pp,ff) = entropy( double(r) );
                        % cov: cross covariance
                        cov_pf = mean2( ( p - p_mean ) .* (f - f_mean ) );
                        c_cov(pp,ff) = - cov_pf;
                        % corr: cross correlation
                        c_corr(pp,ff) = - cov_pf / ( p_std * f_std )
                        % ssim: tructural similarity index
                        c1 = ( 0.01 * L )^2;
                        c2 = ( 0.03 * L )^2;
                        c3 = c2 / 2;
                        c_l(pp,ff) = ( 2 * p_mean * f_mean + c1 ) / ( p_mean^2 + f_mean^2 + c1 );
                        c_c(pp,ff) = ( 2 * p_std * f_std + c2 ) / ( p_std^2 + f_std^2 + c2 );
                        c_s(pp,ff) = ( cov_pf + c3) / ( p_std * f_std + c3 );
                        c_ssim(pp,ff) = - c_l(pp,ff) * c_c(pp,ff) * c_s(pp,ff);
                        % ssim-ml: structural similarity index using Matlab
                        c_ssim_ml(pp,ff) = - ssim( proj_roi(:,:,pp), f ); %'DynamicRange', 'uint16'
                    end
            end
        end
        PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        
        if visualOutput(1)
            h2 = figure('Name', 'Correlation of raw data and flat fields');
            [~, flat_corr_shift_min_pos_x] =  min ( abs( flat_corr_shift_1), [], 2);
            [~, flat_corr_shift_min_pos_y] =  min ( abs( flat_corr_shift_2), [], 2);
                                    
            switch correlation_method
                case 'shift'
                    m = 2; n = 1;
                    
                    subplot(m,n,1);
                    Y = abs(arrayfun(@(x) (flat_corr_shift_2(x,flat_corr_shift_min_pos_y(x))), 1:num_proj_used));
                    plot(Y, '.')
                    axis  tight
                    title(sprintf('minimal absolute vertical shift along rotation axis'))
                    
                    subplot(m,n,2);
                    plot( arrayfun(@(x) (flat_corr_shift_1(x,flat_corr_shift_min_pos_x(x))), 1:num_proj_used) ,'.')
                    axis  tight
                    title(sprintf('minimal absolute horizontal shift (unused)'))
                                   
                case {'diff', 'std', 'entropy', 'ssim', 'ssim-ml', 'cov', 'corr'}
                    m = 2; n = 1;
                    subplot(m,n,1);
                    %f = @(x) normat( min( x, [], 2));                    
                    f = @(x) normat(x(1,:))';
                    Y = [f(c_diff1_l2), f(c_std), f(c_ent), f(c_cov), f(c_corr), f(c_ssim), f(c_ssim_ml)];
                    %Y = [f(c_diff1_l1), f(c_diff1_l2), f(c_diff2_l1), f(c_diff2_l2)];
                    plot( Y, '-' )
                    legend('iso diff L2', 'ratio std dev', 'ratio entropy', 'cov', 'corr', 'ssim', 'ssim-ml' )
                    %legend('anisotropic L1', 'anisotropic L2', 'isotropic L1', 'isotropic L2')
                    axis tight
                    title(sprintf('correlation measures: projection #1'))
                  
                    subplot(m,n,2);
                    f = @(x) normat(x(1,:))';
                    Y = [ f(c_ssim), f(-c_l), f(-c_c), f(-c_s) ];
                    plot( Y, '-' )
                    legend('ssim', 'luminance', 'contrast', 'structure' )
                    axis tight
                    title(sprintf('SSIM and components: projection #1'))
                    
                    % sorted measures: position and values
                    [c_diff1_l2_val,c_diff1_l2_pos] = sort( c_diff1_l2, 2 );
                    [c_std_val,c_std_pos] = sort( c_std, 2 );    
                    [c_ent_val,c_ent_pos] = sort( c_ent, 2 );                        
                    [c_cov_val,c_cov_pos] = sort( c_cov, 2 );
                    [c_corr_val,c_corr_pos] = sort( c_corr, 2 );    
                    [c_ssim_val,c_ssim_pos] = sort( c_ssim, 2 );
                    [c_ssim_ml_val,c_ssim_ml_pos] = sort( c_ssim_ml, 2 );
            end            
            drawnow
        end
        
        PrintVerbose(verbose, '\nFlat- and dark-field correction.')
        t = toc;
        switch correlation_method
            case {'diff', 'std', 'entropy', 'ssim', 'ssim-ml', 'cov', 'corr'}
                switch correlation_method
                    case 'diff'
                        corr_mat = c_diff1_l2;
                    case 'std'
                        corr_mat = c_std;
                    case 'entropy'
                        corr_mat = c_ent;
                    case 'ssim'
                        corr_mat = c_ssim;
                    case 'ssim-ml'
                        corr_mat = c_ssim_ml;
                    case 'cov'
                        corr_mat = c_cov;
                    case 'corr'
                        corr_mat = c_corr;
                end                
                [corr_mat_val, corr_mat_pos] = sort( normat( corr_mat ), 2);
                % Save correlation matrix
                CheckAndMakePath( flatcor_path )
                save( sprintf( '%s/corr_mat_val.mat', flatcor_path), 'corr_mat_val' )
                save( sprintf( '%s/corr_mat_pos.mat', flatcor_path), 'corr_mat_pos' )
                flat_ind = corr_mat_pos(:,1:corr_num_flats);
                % Flat field correction
                parfor nn = 1:num_proj_used                    
                    proj(:, :, nn) = proj(:, :, nn) ./ squeeze( mean( flat(:, :, flat_ind(nn,:)), 3) );
                end
                
            case 'shift'
                % best match
                if corr_shift_max_pixelshift == 0
                    corr_mat = flat_corr_shift_2;
                    [~, pos] = min( abs( corr_mat ), [], 2 );
                    parfor nn = 1:num_proj_used
                        proj(:, :, nn) = proj(:, :, nn) ./ flat(:, :, pos(nn));
                    end
                    
                    % use all flats which are shifted less pixels than corr_shift_max_pixelshift
                elseif corr_shift_max_pixelshift > 0
                    nflats = zeros(1, num_proj_used);
                    parfor nn = 1:num_proj_used
                        vec = 1:num_ref_used;
                        flat_ind = vec( abs( flat_corr_shift_2(nn, :) ) < corr_shift_max_pixelshift );
                        if numel( flat_ind ) > corr_num_flats
                            flat_ind( corr_num_flats + 1:end ) = [];
                        end
                        if isempty( flat_ind )
                            [~, flat_ind] = min( abs( flat_corr_shift_2(nn, :) ) );
                        end
                        nflats(nn) = numel(flat_ind);
                        proj(:, :, nn) = proj(:, :, nn) ./ squeeze( mean( flat(:, :, flat_ind), 3) );
                    end
                else
                    error('Value of maximum shift (%g) is not >= 0', corr_shift_max_pixelshift)
                end
        end
        
        PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        if strcmp( correlation_method, 'shift') && corr_shift_max_pixelshift > 0
            PrintVerbose(verbose, '\n number of flats used per projection: [mean, min, max] = [%g, %g, %g]', mean( nflats ), min( nflats ), max( nflats) )
        end
    end
    PrintVerbose(verbose, '\n sinogram size = [%g, %g, %g]', size( proj ) )
    if visualOutput(1)
        figure(h1)
        subplot(2,2,4)
        imsc1( proj(:,:,1))
        axis equal tight
        title(sprintf('flat-&-dark corrected projection #1'))
        colorbar
        drawnow
    end
    
    %% Ring artifact filter
    sino_slice = round( size( proj, 2) / 2 );    
    sino_unfilt = squeeze( proj(:,sino_slice,:) )';
    if ring_filter(1)
        t = toc;
        PrintVerbose(verbose, '\nFilter ring artifacts.')
        switch lower( ring_filter_method )
            
            % Combined wavelet-FFT filter to remove ring artrifacts
            case 'wavelet-fft'                
                
                parfor nn = 1:size( proj, 2)
                    sino = squeeze( proj(:,nn,:) )';
                    [d2, d1] = size( sino );
                    sino = FilterStripesCombinedWaveletFFT( sino, dec_levels, wname, sigma )';
                    sino = sino( 1:d1, 1:d2 );
                    proj(:,nn,:) = sino;
                end
                PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
               
                if visualOutput(1)
                    
                    h3 = figure('Name', 'Sinogram and ring filter');
                    
                    subplot(3,1,1)
                    [~, sorted_angle_index] = sort( angles );
                    imsc( sino_unfilt(:,sorted_angle_index) )
                    axis equal tight
                    title(sprintf('sino unfiltered, y = %u', sino_slice))
                    colorbar 

                    subplot(3,1,2)
                    imsc( squeeze( proj(:,sino_slice,:) )' )
                    axis equal tight
                    title(sprintf('sino filtered, y = %u', sino_slice))
                    colorbar
                    
                    subplot(3,1,3)
                    imsc( ( squeeze( proj(:,sino_slice,:) )' - sino_unfilt ) )
                    axis equal tight
                    title(sprintf('sino filt - sino, y = %u', sino_slice))
                    colorbar
                    
                    drawnow
                end
                
            % Simple ring artifact removal filter
            case 'jm'                                                
                if numel( ring_filter_median_width ) > 1
                    %% Combine if/else
                    for nn = ring_filter_median_width
                        %m=sum(im);
                        proj_mean = sum( proj, 3);
                        
                        %md=medfilt1(m,nn);
                        proj_mean_med = medfilt2( proj_mean, [nn, 1], 'symmetric' );
                                                
                        %f=md./m;
                        mask = proj_mean_med ./ proj_mean;
                                                
                        %im=bsxfun(@times,im,f);
                        proj = bsxfun( @times, proj, mask);                                                
                    end                    
                else                    
                    proj_mean = mean( proj, 3);
                    proj_mean_med = medfilt2( proj_mean, [ring_filter_median_width, 1], 'symmetric' );
                    mask = proj_mean_med ./ proj_mean;
                    proj = bsxfun( @times, proj, mask);
                end
                
                PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
                PrintVerbose( verbose, '\n ring filter mask min/max: %f, %f', min( mask(:) ), max( mask(:) ) )
                if visualOutput(1)
                    h3 = figure('Name', 'Sinogram and ring filter');
                                        
                    subplot(2,2,1)
                    imsc( sino_unfilt )
                    axis equal tight
                    title(sprintf('sino unfiltered, y = %u', sino_slice))
                    colorbar %('Location', 'southoutside')
                    
                    subplot(2,2,2)
                    imsc( ( squeeze( proj(:,sino_slice,:) ) - sino_unfilt' )' )
                    axis equal tight
                    title(sprintf('sino filt - sino, y = %u', sino_slice))
                    colorbar %('Location', 'southoutside')
                    
                    subplot(2,2,3)
                    imsc1( proj_mean )
                    axis equal tight
                    title(sprintf('mean projection'))
                    colorbar %('Location', 'southoutside')
                    
                    subplot(2,2,4)
                    imsc1( FilterHisto( mask, 5 ) )
                    axis equal tight
                    title(sprintf('mask for normalization'))
                    colorbar %('Location', 'southoutside')
                    
                    drawnow
                end
        end
    end
    proj_min = min( proj(:) );
    proj_max = max( proj(:) );
    
    %% Write corrected projections
    if write_flatcor(1)
        t = toc;
        PrintVerbose(verbose, '\nSave flat-corrected projections.')
        CheckAndMakePath( flatcor_path )
        % Projections
        parfor nn = 1:num_proj_used
            filename = sprintf('%sproj_%06u.tif', flatcor_path, nn );
            write32bitTIFfromSingle(filename, rot90( proj(:, :, nn) ) );
        end
        PrintVerbose(verbose, ' Time elapsed: %.1f (%.2f min)', toc-t, (toc-t)/60)
    end
end

%% Rotation axis position and tomgraphic reconstruction parameters %%%%%%%%
if do_phase_retrieval(1)
    if isempty( take_neg_log )
        take_neg_log = 0;
    end
end
tint = 0;
if do_tomo(1)
    t = toc;
    
    % ROI for correlation of projections at angles 0 & pi
    if isempty( rot_corr_area1 )
        switch excentric_rot_axis
            case -1
                rot_corr_area1 = [0 0.25];
            case 0
                rot_corr_area1 = [0.25 0.75];
            case 1
                rot_corr_area1 = [0.75 1];
        end
    end
    if isempty( rot_corr_area2 )
        rot_corr_area2 = [0.1 0.9];
    end
    rot_corr_area1 = IndexParameterToRange( rot_corr_area1, raw_im_shape_binned1 );
    rot_corr_area2 = IndexParameterToRange( rot_corr_area2, raw_im_shape_binned2 );
    
    % Full rotation angle
    if isempty( rot_angle_full )
        if isfield( par, 'rotation')
            % From log file
            rot_angle_full = par.rotation / 180 * pi;
        elseif exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
            % from beam current log
            rot_angle_full = (cur.proj(end).angle - cur.proj(1).angle) * pi /180; % KIT: , EHD: ok
        else
            % Guess from correlation of projections
            im1 = proj( rot_corr_area1, rot_corr_area2, 1);
            im2 = proj( rot_corr_area1, rot_corr_area2, num_proj_used);
            [~, cm1] = ImageCorrelation( im1, im2, 0, 0, 0); % large if full angle is  2 * pi
            [~, cm2] = ImageCorrelation( im1, flipud(im2), 0, 0, 0); % large if full angle is pi
            if max( cm1(:) ) > max( cm2(:) )
                rot_angle_full = 2 * pi;
            elseif max( cm1(:) ) < max( cm2(:) )
                rot_angle_full = pi;
            else
                error('Determination of full angle of rotation from correlation of first and last projection not successful.')
            end
        end
    end
    PrintVerbose(verbose, '\n full rotation angle: %g * pi', rot_angle_full / pi)
    if exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
        % for KIT cam this includes missing angles
        angles = [cur.proj.angle] / 180 * pi;
        if strcmpi(cam, 'kit')
            % drop angles where projections are missing
            angles = angles(1 + proj_nums);
        else
            angles = angles(proj_range);
        end
    else
        if isfield( par, 'projections' )
            num_proj = double( par.projections );
        elseif isfield( par, 'n_angles' )
            num_proj = double( par.n_angles );
        end
        switch lower( cam )
            case 'ehd'
                angles = rot_angle_full * (0:num_proj - 1) / (num_proj - 1); % EHD: ok
            case 'kit'
                angles = rot_angle_full * (0:num_proj - 1) / num_proj; % KIT: ok if par.projections exist
        end
    end
    if numel( angles ) ~= num_proj_used
        error('Number of elements in array of angles (%g) unequal number of projections read (%g)', numel( angles ), num_proj_used)
    end
    % retrieve index at angles 0 and pi
    [val1, ind1] = min( abs( angles ));
    [val2, ind2] = min( abs( angles - pi ));
    % ROI
    im1 = proj( rot_corr_area1, rot_corr_area2, ind1);
    im2 = flipud( proj( rot_corr_area1, rot_corr_area2, ind2) );
    PrintVerbose(verbose, '\n correlation of images : [filename index, projection index, angle / pi] = [%g, %g, %g] and [%g, %g, %g]', ... 
        proj_nums(ind1), ind1, val1 / pi, proj_nums(ind2), ind2, (val2 + pi) / pi )
    if rot_corr_gradient(1)
        l = 2;
        [g11, g12] = gradient(im1);
        im1 = abs(g11).^l + abs(g12).^l;
        [g21, g22] = gradient(im2);
        im2 = abs(g21).^l + abs(g22).^l;
    end
    % Correlation
    out = ImageCorrelation( im1, im2, 0, 0, 0);
    % relative shift
    rot_corr_shift_x = round( out.shift1, decimal_round_precision) + rot_corr_area1(1) + (rot_corr_area1(end) - raw_im_shape_binned1) - 1;
    rot_corr_shift_y = round( out.shift2, decimal_round_precision) + rot_corr_area2(1) + (rot_corr_area1(end) - raw_im_shape_binned1) - 1;
    PrintVerbose(verbose, '\n relative shift: %g, %g', rot_corr_shift_x, rot_corr_shift_y)
    rot_axis_offset_calc = rot_corr_shift_x / 2;
    rot_axis_pos_calc = raw_im_shape_binned1 / 2 + rot_axis_offset_calc;
    PrintVerbose(verbose, '\n calculated rotation axis offset: %.2f', rot_axis_offset_calc)
    % use calculated offset if both offset and position are empty
    if isempty( rot_axis_offset ) && isempty( rot_axis_pos )
        rot_axis_offset = rot_axis_offset_calc;
    end
    if isempty( rot_axis_pos )
        rot_axis_pos = raw_im_shape_binned1 / 2 + rot_axis_offset;
    end
    if isempty( rot_axis_offset )
        rot_axis_offset = rot_axis_pos - raw_im_shape_binned1 / 2;
    end
    
    % Tilt of rotation axis
    if interactive_determination_of_rot_axis_tilt(1)
        im1c = RotAxisSymmetricCropping( proj(:,rot_corr_area2,ind1), rot_axis_pos, 1);
        im2c = flipud(RotAxisSymmetricCropping( proj(:,rot_corr_area2,ind2) , rot_axis_pos, 1));
        [optimizer, metric] = imregconfig('monomodal');
        tform = imregtform(im2c, im1c, 'rigid', optimizer, metric);
        rot_axis_tilt_calc = asin( tform.T(1,2) ) / 2;
        PrintVerbose(verbose, '\n calculated tilt of rotation axis : %g rad = %g deg', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
    else
        rot_axis_tilt_calc = [];
    end
    if isempty( rot_axis_tilt )
        if abs( rot_axis_pos_calc ) < 0.01 % 0.573 degree
            rot_axis_tilt = rot_axis_tilt_calc;
            PrintVerbose(verbose, '\n use calculated tilt of rotation axis')
        else
            rot_axis_tilt = 0;
        end
    end
    
    %% Determine rotation axis position
    tint = 0;
    if interactive_determination_of_rot_axis(1)
        tint = toc;        
        fprintf( '\n\nENTER INTERACTIVE MODE' )
        fprintf( '\n number of pixels: %u', raw_im_shape_binned1)
        fprintf( '\n image center: %.1f', raw_im_shape_binned1 / 2)
        
        fra_take_neg_log = 1;
        fra_number_of_stds = 4;
        fra_vol_shape = [];
        
        if interactive_determination_of_rot_axis_slice > 1
            slice = interactive_determination_of_rot_axis_slice;
        elseif interactive_determination_of_rot_axis_slice <= 1 && interactive_determination_of_rot_axis_slice >= 0
            slice = round((raw_im_shape_binned2 - 1) * interactive_determination_of_rot_axis_slice + 1 );
        end
        fprintf( '\n slice : %u', slice)
        
        fprintf( '\n\nOFFSET:' )
        fprintf( '\n current rotation axis offset / position : %.2f, %.2f', rot_axis_offset, rot_axis_pos)
        fprintf( '\n calcul. rotation axis offset / position : %.2f, %.2f', rot_axis_offset_calc, rot_axis_pos_calc)
        fprintf( '\n default offset range : current ROT_AXIS_OFFSET + (-4:0.5:4)')
        offset = input( '\n\nENTER RANGE OF ROTATION AXIS OFFSETS (if empty use default range, if scalar skips interactive mode): ');
        if isempty( offset )
            % default range is centered at the given or calculated offset
            offset = rot_axis_offset + (-4:0.5:4);
        end                   
        
        % Loop over offsets
        while ~isscalar( offset )
            
            % Reco
            [vol, metrics_offset] = find_rot_axis_offset(proj, angles, slice, offset, rot_axis_tilt, fra_take_neg_log, fra_number_of_stds, fra_vol_shape);
            
            % Metric minima
            [~, min_pos] = min(cell2mat({metrics_offset(:).val}));
            [~, max_pos] = max(cell2mat({metrics_offset(:).val}));
            
            % Print image number, rotation axis values, and different metrics
            fprintf( '\n\nOFFSET:' )
            fprintf( '\n current rotation axis offset/position : %.2f, %.2f\n', rot_axis_offset, rot_axis_pos)
            fprintf( '%11s', 'image no.', 'offset', metrics_offset.name)
            for nn = 1:numel(offset)
                if offset(nn) == rot_axis_offset
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
            x = [1:4 6:7];
            Y = cell2mat({metrics_offset(x).val});
            plot( offset, Y, '-+');
            axis tight
            legend( metrics_offset(x).name)            
            title(sprintf('metric VS rotation axis offset'))
            drawnow            
            
            % Play
            nimplay(vol, 0, 0, 'OFFSET: sequence of reconstructed slices using different rotation axis offsets')
            
            % Input
            offset = input( '\n\nENTER ROTATION AXIS OFFSET OR A RANGE OF OFFSETS (if empty use current offset): ');
            if isempty( offset )
                offset = rot_axis_offset;
            end
            
            if isscalar( offset )
                fprintf( ' old rotation axis offset : %.2f', rot_axis_offset)
                rot_axis_offset = offset;
                fprintf( '\n new rotation axis offset : %.2f', rot_axis_offset)

               if interactive_determination_of_rot_axis_tilt(1)
                
                % Tilt
                fprintf( '\n\nTILT:' ) 
                fprintf( '\n current rotation axis tilt : %g rad = %g deg', rot_axis_tilt, rot_axis_tilt * 180 / pi)
                fprintf( '\n calcul. rotation axis tilt : %g rad = %g deg', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
                fprintf( '\n default tilt range is : current ROT_AXIS_TILT + (-0.005:0.001:0.005)')
                tilt = input( '\n\nENTER TILT OF ROTATION AXIS OR RANGE OF TILTS (if empty use default):');        
                if isempty( tilt )
                    % default range is centered at the given or calculated tilt
                    tilt = rot_axis_tilt + (-0.005:0.001:0.005);                    
                end
                
                while ~isscalar( tilt )
                    
                    % Reco
                    [vol, metrics_tilt] = find_rot_axis_tilt( proj, angles, slice, offset, tilt, fra_take_neg_log, fra_number_of_stds, fra_vol_shape);
                    
                    % Metric minima
                    [~, min_pos] = min(cell2mat({metrics_tilt(:).val}));
                    [~, max_pos] = max(cell2mat({metrics_tilt(:).val}));
                    
                    % Print image number and rotation axis tilt
                    fprintf( '%11s', 'image no.', 'tilt/rad', 'tilt/deg', metrics_tilt.name )
                    for nn = 1:numel(tilt)
                        if tilt(nn) == rot_axis_tilt
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
                    legend( metrics_tilt(x).name)
                    title(sprintf('metric VS rotation axis tilt'))
                    drawnow
                    
                    
                    % Play
                    nimplay(vol, 0, 0, 'TILT: sequence of reconstructed slices using different rotation axis tilts')
                    
                    % Input
                    tilt = input( '\nENTER TILT OF ROTATION AXIS OR RANGE OF TILTS (if empty use current tilt):');
                    if isempty( tilt )
                        tilt = rot_axis_tilt;
                    end
                    
                    if isscalar( tilt )
                        rot_axis_tilt = tilt;
                                                                                               
                        rot_axis_pos = raw_im_shape_binned1 / 2 + rot_axis_offset;
                        
                        % Compare projection at 0 pi and projection at 1 pi corrected for rotation axis tilt                        
                        im1c = RotAxisSymmetricCropping( proj(:,rot_corr_area2,ind1), rot_axis_pos, 1);
                        im2c = flipud(RotAxisSymmetricCropping( proj(:,rot_corr_area2,ind2) , rot_axis_pos, 1));
                        [optimizer, metric] = imregconfig('monomodal');
                        tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);
                        rot_axis_tilt_calc = asin( tform_calc.T(1,2) ) / 2;
                        im2c_warped_calc =  imwarp(im2c, tform_calc, 'OutputView', imref2d(size(im1c)));
                              
                        tform_cur = tform_calc;
                        %tform_cur.Dimensionality = 2;
                        tform_cur.T = [cos( 2 * rot_axis_tilt ) sin( 2 * rot_axis_tilt ) 0; ...
                                      -sin( 2 * rot_axis_tilt ) cos( 2 * rot_axis_tilt ) 0 ; ...
                                      tform_calc.T(3,1) tform_calc.T(3,2) 1];
                                      %0 0 1];
                                      
                        im2c_warped_cur =  imwarp(im2c, tform_cur, 'OutputView', imref2d(size(im1c)));
                        
                        x = ceil( 2 * abs( sin(rot_axis_tilt) ) * max( size(im1c)) ) + 2;
                        
                        im2c_warped_cur(im2c_warped_cur == 0) = mean2( im2c );
                        im2c_warped_calc(im2c_warped_calc == 0) = mean2( im2c );                        
                        
                        fprintf( '\n current rotation axis tilt: %g rad (%g deg)', rot_axis_tilt, rot_axis_tilt * 180 / pi)
                        fprintf( '\n calcul. rotation axis tilt: %g rad (%g deg)', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)                        
                                                                        
                        name = sprintf( 'projections at 0 and 180 degree. corrected. CURRENT rot axis tilt: %g, rot axis offset: %g', rot_axis_tilt, rot_axis_offset);
                        nimplay( cat(3, im1c(x:end-x,x:end-x)', im2c_warped_cur(x:end-x,x:end-x)'), 0, 0, name)
                                                
                        name = sprintf( 'projections at 0 and 180 degree. corrected. CALCULATED rot axis tilt: %g, rot axis offset: %g', rot_axis_tilt_calc, rot_axis_offset);
                        nimplay( cat(3, im1c(x:end-x,x:end-x)', im2c_warped_calc(x:end-x,x:end-x)'), 0, 0, name)
                                                                                                
                        tilt = input( '\nENTER ROTATION AXIS TILT (if empty use current value): ');
                        if isempty( tilt )
                            tilt = rot_axis_tilt;
                        else
                            rot_axis_tilt = tilt;
                        end
                        
                        offset = input( '\nENTER RANGE OF OFFSETS TO CONTINUE INTERACTIVE LOOP OR TYPE ENTER TO EXIT LOOP: ');
                        if isempty( offset )
                            offset = rot_axis_offset;
                        end                    
                    end                    
                end
                else
                    rot_axis_tilt = 0;
                end
            end
        end
             
        rot_axis_pos = raw_im_shape_binned1 / 2 + rot_axis_offset;
        fprintf( '\nEND OF INTERACTIVE MODE\n' )
   
        tint = toc - tint;
    end   
    PrintVerbose(verbose, '\n rotation axis offset: %.2f', rot_axis_offset );
    PrintVerbose(verbose, '\n rotation axis position: %.2f', rot_axis_pos );
    PrintVerbose(verbose, '\n rotation axis tilt: %g rad (%g deg)', rot_axis_tilt, rot_axis_tilt * 180 / pi)
       
    if interactive_determination_of_rot_axis_tilt(1) && visualOutput(1)
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
        
        x = ceil(1 + 2 * abs( sin(rot_axis_tilt) ) * max( size(im1c)) );
        
        subplot(m, n, 3)
        imsc1( abs( im1c(x:end-x,x:end-x) - im2c(x:end-x,x:end-x) ) )
        axis equal tight
        title(sprintf('difference before'))
        colorbar
        
        subplot(m, n, 4)
        
        %% tform_cal or tform_cur
        [optimizer, metric] = imregconfig('monomodal');
        tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);
        im2c_warped_cur =  imwarp(im2c, tform_calc, 'OutputView',imref2d(size(im1c)));
        imsc1( abs( im1c(x:end-x,x:end-x) - im2c_warped_cur(x:end-x,x:end-x) ) )
        axis equal tight
        title(sprintf('difference corrected'))
        colorbar
        
        drawnow              
    end
end

%% Stitch projections
if stitch_projections(1)
    t = toc;
    PrintVerbose(verbose, '\nStitch projections:')
    if rot_angle_full < 1.9 * pi
        error( 'full angle of rotation smaller than 2 pi: %g pi', rot_angle_full/pi)
    end
    % number of unstitched projections
    num_proj = size( proj, 3);
    % last projection within [0,pi)
    [~, num_proj_sti] = min( abs(angles - pi));
    % number of stitched projections
    num_proj_sti = num_proj_sti - 1;
    % index range of projections to be stitched
    xl = 1:round(rot_axis_pos);
    xr = 1:xl(end)-1;
    im_shape_sti1 = numel( xl ) + numel( xr );
    % Preallocation
    proj_sti = zeros( im_shape_sti1 , raw_im_shape_binned2, num_proj_sti, 'single');
    for nn = 1:num_proj_sti
        nn2 = mod(num_proj_sti + nn - 1, num_proj) + 1;
        im = zeros( im_shape_sti1, raw_im_shape_binned2);
        switch lower( stitch_method )
            case 'step'
                im = cat(1, proj(xl,:,nn), flipud( proj(xr,:,nn2) ) );
            case {'linear', 'sine'}
                % overlap region
                overlap = 2 * rot_axis_pos - raw_im_shape_binned1 : raw_im_shape_binned1;
                % overlap ramp
                x = (0:1/(numel(overlap)-1):1);
                % 1D weight
                w = ones(raw_im_shape_binned1, 1);
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
                im(1:raw_im_shape_binned1,:) = iml;
                im(end - raw_im_shape_binned1 + 1:end,:) = im(end - raw_im_shape_binned1 + 1:end,:) + imr;
        end
        proj_sti(:,:,nn) = im;
    end
    pause(0.01)
    proj = proj_sti;
    clear proj_sti;
    angles = angles(1:num_proj_sti);
    PrintVerbose(verbose, ' Time elapsed: %.1f (%.2f min)', toc-t, (toc-t)/60)
    PrintVerbose(verbose, '\n shape of stitched projections : %u %u %u', size( proj ) )
    PrintVerbose(verbose, '\n memory allocted : %.2 GiB', Bytes( proj, 3 ) )
end

%% Crop projections at rotation axis position
% to avoid oversampling for scans with excentric rotation axis and
% reconstructing without stitching
if crop_at_rot_axis(1)
    switch excentric_rot_axis
        case 1
            proj( ceil(rot_axis_pos) + 1:end, :, :) = [];
        case -1
            %% CHECK
            proj( 1:floor(rot_axis_pos)-1, :, :) = [];
    end
    if isempty( vol_shape )                
        vol_shape = [raw_im_shape_binned1, raw_im_shape_binned1, raw_im_shape_binned2];        
    end
end

%% Save sinograms
if write_sino(1)
    t = toc;
    PrintVerbose(verbose, '\nSave sinogram:')
    CheckAndMakePath(sino_path)
    % Save slices
    parfor nn = 1:raw_im_shape_binned2
        filename = sprintf( '%ssino_%06u.tif', sino_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    pause(0.01)
    PrintVerbose(verbose, ' done.')
    PrintVerbose(verbose, ' Time elapsed: %g s (%.2f min)', toc-t, (toc-t)/60)
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edp = [energy, sample_detector_distance, eff_pixel_size_binned];
if do_phase_retrieval(1)
    PrintVerbose(verbose, '\nPhase retrieval.')
    t = toc;
    PrintVerbose(verbose, '\n energy : %g eV', energy)
    PrintVerbose(verbose, '\n sample detector distance : %g m', sample_detector_distance)
    
    if isempty( take_neg_log )
        take_neg_log = 0;
    end
    
    % Phase retrieval filter
    im_shape = [size(proj,1) , size(proj,2)];
    im_shape_pad = (1 + phase_padding) * im_shape;
    [phase_filter, pha_appendix] = PhaseFilter( phase_retrieval_method, im_shape_pad, edp, phase_retrieval_reg_par, phase_retrieval_bin_filt, phase_retrieval_cutoff_frequ, 'single');
    
    % reco phase dir
    if isempty( subfolder_reco )
        reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep];
    else
        reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep, subfolder_reco, filesep];
    end
    CheckAndMakePath( reco_phase_path )
    PrintVerbose(verbose, '\n reco_phase_path : %s', reco_phase_path)
    PrintVerbose(verbose, '\n phase retrieval method : %s', phase_retrieval_method)
    
    % Retrieval
    parfor nn = 1:size(proj, 3)
        % combined GPU and parfor usage requires memory management
        im = padarray( proj(:,:,nn), phase_padding * im_shape, 'symmetric', 'post' );
        %im = padarray( gpuArray( proj(:,:,nn) ), raw_im_shape_binned, 'post', 'symmetric' );
        im = -real( ifft2( phase_filter .* fft2( im ) ) );
        proj(:,:,nn) = im(1:im_shape(1), 1:im_shape(2));
        %proj(:,:,nn) = gather( im(1:raw_im_shape_binned1, 1:raw_im_shape_binned2) );
    end
    pause(0.01)
    PrintVerbose(verbose, ' Time elapsed: %g s (%.2f min)', toc-t, (toc-t)/60)
    
    % Save phase maps
    if write_phase_map(1)
        t = toc;
        PrintVerbose(verbose, '\nSave phase maps:')
        phase_map_path = [phase_map_path, pha_appendix, filesep];
        CheckAndMakePath( phase_map_path );
        parfor nn = 1:size( proj, 3)
            filename = sprintf( '%sphase_%06u.tif', phase_map_path, nn);
            write32bitTIFfromSingle( filename, squeeze( rot90( proj( :, :, nn) ) ) )
        end
        pause(0.01)
        PrintVerbose(verbose, ' Time elapsed: %g s (%.2f min)', toc-t, (toc-t)/60)
    end
end

%% Save sinograms of phase maps
if write_sino_phase(1)
    t = toc;
    PrintVerbose(verbose, '\nSave phase map sinogram:')
    CheckAndMakePath(sino_phase_path)
    % Save slices
    parfor nn = 1:raw_im_shape_binned2
        filename = sprintf( '%ssino_%06u.tif', sino_phase_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    pause(0.01)
    PrintVerbose(verbose, ' done.')
    PrintVerbose(verbose, ' Time elapsed: %g s (%.2f min)', toc-t, (toc-t)/60)
end

%% Tomographic reco
if do_tomo(1)
    
    if stitch_projections(1)
        rot_axis_offset_reco = 0;
    elseif crop_at_rot_axis(1)
        rot_axis_offset_reco = rot_axis_pos - size( proj, 1) / 2;
    else
        rot_axis_offset_reco = rot_axis_offset;
    end
    
    %% Volume shape and size
    if isempty( vol_shape )
        % default volume given by the detector width and height
        vol_shape1 = size( proj, 1);
        vol_shape = [vol_shape1, vol_shape1, raw_im_shape_binned2];
    else
        vol_shape1 = size( proj, 1);
        if vol_shape(1) <=  1
            vol_shape(1) = round( 1 + vol_shape(1) * (vol_shape1 - 1) );
        end
        if vol_shape(2) <=  1
            vol_shape(2) = round( 1 + vol_shape(2) * (vol_shape1 - 1) );
        end
        if vol_shape(3) <=  1
            vol_shape(3) = round( 1 + vol_shape(3) * (raw_im_shape_binned2 - 1) );
        end
    end
    if isempty( vol_size )
        vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];
    end
    PrintVerbose(verbose, '\n shape reconstructed volume: [%g, %g, %g]', vol_shape )    
    
    PrintVerbose(verbose, '\nTomographic reconstruction:')
    if isempty( take_neg_log )
        take_neg_log = 1;
    end
    vol_min = NaN;
    vol_max = NaN;
    angles_reco = angles;
    % Delete redundant projection
    if isequal( angles(1), angles(end) )
        angles_reco(end) = [];
        proj(:,:,end) = [];
    end
    
    % Filter sinogram
    PrintVerbose(verbose, '\n Filter sino:' )
    t2 = toc;
    filt = iradonDesignFilter(fbp_filter_type, (1 + fbp_filter_padding) * size( proj, 1), fpb_filter_freq_cutoff);
    if butterworth_filter(1)
        [b, a] = butter(butterworth_order, butterworth_cutoff_frequ);
        bw = freqz(b, a, numel(filt) );
        filt = filt .* bw;
    end
    
    proj_shape1 = size( proj, 1);
    parfor nn =  1:size( proj, 2)
        im = proj(:,nn,:);
        im = padarray( NegLog(im, take_neg_log), fbp_filter_padding * [proj_shape1 0 0], fbp_filter_padding_method, 'post' );        
        im = real( ifft( bsxfun(@times, fft( im, [], 1), filt), [], 1, 'symmetric') ); %% TODO: Check symmetric option
        %im = real( ifft( bsxfun(@times, fft( im, [], 1), filt), [], 1) );        
        proj(:,nn,:) = im(1:proj_shape1,:,:);
    end    
    pause(0.01)
    PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
      
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
    PrintVerbose(verbose, '\n Backproject:')
    t2 = toc;
    vol = astra_parallel3D( permute(proj, [1 3 2]), rot_angle_offset + angles_reco, rot_axis_offset_reco, vol_shape, vol_size, astra_pixel_size, link_data, rot_axis_tilt, gpu_index);
    pause(0.01)    
    PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
    
    % Save volume
    if write_reco(1)        
        if do_phase_retrieval(1)
            reco_path = reco_phase_path;
        end
        CheckAndMakePath(reco_path)
        
        % Save reco path and scan path in file
        filename = [userpath, filesep, 'experiments/p05/pathtolastreco'];
        fid = fopen( filename , 'w' );
        fprintf( fid, '%s', reco_path );
        fclose( fid );
        
        % Single precision: 32-bit float tiff        
        if write_float(1)
            PrintVerbose(verbose, '\n Write floats:')
            t2 = toc;
            save_path = sprintf( '%sfloat_rawBin%u/', reco_path, raw_bin);
            CheckAndMakePath(save_path)
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_path, nn);
                write32bitTIFfromSingle( filename, vol( :, :, nn) )
            end
            pause(0.01)
            PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
        end
        
        % Binned single precision: 32-bit float tiff
        if write_float_binned(1)
            PrintVerbose(verbose, '\n Write float binned:')
            t2 = toc;
            save_path = sprintf( '%sfloat_rawBin%u_recoBin%u/', reco_path, raw_bin, reco_bin);
            CheckAndMakePath(save_path)
            for ii = 1:floor( size( vol, 3) / reco_bin )
                filename = sprintf( '%sreco_%06u.tif', save_path, ii);
                nn = 1 + reco_bin * (ii - 1) + (0:reco_bin - 1);                
                im = Binning( sum(vol( :, :, nn), 3 ), reco_bin) / reco_bin^3;
                write32bitTIFfromSingle( filename, im )
            end
            pause(0.01)
            PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
        end
        
        vol_min = min( vol(:) );
        vol_max = max( vol(:) );
        
        % Compression of dynamic range
        if write_8bit || write_16bit || write_8bit_binned
            
            % Compression method
            switch compression
                case 'full'
                    % use the full dynamic range of the data
                    vol = (vol - vol_min) ./ (vol_max - vol_min);
                    tlow = vol_min;
                    thigh = vol_max;
                case 'std'
                    % use the dynamic range within several standard
                    % deviations centered around the mean value of the data
                    vol_mean = mean( vol(:) );
                    vol_std = std( vol(:) );
                    tlow = vol_mean - compression_std_num * vol_std;
                    thigh = vol_mean + compression_std_num * vol_std;                    
                    % Rescale volume
                    vol = (vol - tlow) / (thigh - tlow);
                case 'threshold'
                    % use dyanmic range within given lower and upper
                    % thresholds
                    tlow = compression_threshold(1);
                    thigh = compression_threshold(2);
                    vol = ( vol - tlow ) / ( thigh - tlow );
                case 'histo'
                    % crop dynamic range using histogram
                    t2 = toc;
                    PrintVerbose(verbose, '\nCompression by histogram thresholds:' )
                    num_bins_max = 2048;
                    bl = min( round( 1 / compression_histo(1) / 2 ) * 2, num_bins_max);
                    bh = min( round( 1 / compression_histo(2) / 2 ) * 2, num_bins_max);
                    num_bins = min( lcm( bl, bh ), num_bins_max);
                    PrintVerbose(verbose, '\n lower/higher percentage : %g%%, %g%%', compression_histo*100 )
                    PrintVerbose(verbose, '\n bins required for lower/higher percentage thresholds : %u, %u', bl, bh )
                    PrintVerbose(verbose, '\n number of histogram bins (least common multiple) : %u', num_bins )
                    [x, y, z] = size( vol );
                    xx = round( x / 2 ) + (-ceil(0.9*x/2/sqrt(2)):floor(0.9*x/2/sqrt(2)));
                    yy = round( y / 2 ) + (-ceil(0.9*y/2/sqrt(2)):floor(0.9*y/2/sqrt(2)));
                    zz = round( z / 2 ) + (-ceil(z/4):floor(z/4));
                    vol_min = min3( vol(xx,yy,zz) );
                    vol_max = max3( vol(xx,yy,zz) );
                    PrintVerbose(verbose, '\n volume ROI  min/max : %g, %g', vol_min, vol_max )
                    [histo, edges] = histcounts( vol(xx,xx,zz), num_bins, 'BinLimits', [vol_min vol_max] );
                    tlow = edges( round( num_bins / bl ) );
                    thigh = edges( end - round( num_bins / bh ) );
                    PrintVerbose(verbose, ' \n lower / higher threshold : %g, %g', tlow, thigh )
                    % Normalize volume
                    vol = (vol - tlow) / (thigh - tlow);                                   
                    PrintVerbose(verbose, '\n done in %.0f s (%.2f min).',toc - t2, (toc - t2) / 60)
            end
            
            % 16-bit tiff
            if write_16bit(1)
                PrintVerbose(verbose, '\n Write uint16:')
                t2 = toc;
                save_path = sprintf( '%s/uint16_rawBin%u/', reco_path, raw_bin);
                CheckAndMakePath(save_path)
                parfor ii = 1:size( vol, 3)
                    filename = sprintf( '%sreco_%06u.tif', save_path, ii);
                    imwrite( uint16( (2^16 - 1) * vol( :, :, ii) ), filename );
                end
                pause(0.01)                
                PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
            end
            
            % 8-bit tiff
            if write_8bit(1)
                PrintVerbose(verbose, '\n Write uint8:')
                t2 = toc;
                save_path = sprintf( '%s/uint8_rawBin%u/', reco_path, raw_bin);                
                CheckAndMakePath(save_path)
                parfor ii = 1:size( vol, 3)
                    filename = sprintf( '%sreco_%06u.tif', save_path, ii);
                    imwrite( uint8( (2^8 - 1) * vol( :, :, ii) ), filename );
                end
                pause(0.01)                
                PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
            end
            
            % 8-bit tiff binned
            if write_8bit_binned(1)
                PrintVerbose(verbose, '\n Write uint8 binned:')
                t2 = toc;
                save_path = sprintf( '%s/uint8_rawBin%u_recoBin%u/', reco_path, raw_bin, reco_bin);
                CheckAndMakePath(save_path)
                for ii = 1:floor( size( vol, 3) / reco_bin )
                    filename = sprintf( '%sreco_%06u.tif', save_path, ii);
                    nn = 1 + reco_bin * (ii - 1) + (0:reco_bin - 1);
                    im =  Binning( sum(vol( :, :, nn), 3 ), reco_bin) / reco_bin^3;
                    imwrite( uint8( (2^8 - 1) * im ), filename );
                end
                pause(0.01)                
                PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
            end          
        end
    end
    
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60 )
end
    
%% Log file
if write_reco(1)
    logfile_path = reco_path;
else
    logfile_path = out_path;
end
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
fprintf(fid, 'raw_image_shape_raw : %u %u\n', raw_im_shape_raw);
fprintf(fid, 'raw_roi : %u %u\n', raw_roi);
fprintf(fid, 'raw_image_shape : %u %u\n', raw_im_shape);
fprintf(fid, 'raw_image_shape_binned : %u %u\n', raw_im_shape_binned);
fprintf(fid, 'raw_binning_factor : %u\n', raw_bin);
fprintf(fid, 'effective_pixel_size : %g micron\n', eff_pixel_size * 1e6);
fprintf(fid, 'effective_pixel_size_binned : %g micron\n', eff_pixel_size_binned * 1e6);
fprintf(fid, 'energy : %g eV\n', energy);
fprintf(fid, 'sample_detector_distance : %f m\n', sample_detector_distance);
fprintf(fid, 'norm_by_ring_current : %u\n', norm_by_ring_current);
fprintf(fid, 'proj_flat_correlation_method : %s\n', correlation_method);
fprintf(fid, 'proj_flat_correlation_num_flats : %u\n', corr_num_flats);
fprintf(fid, 'flat_field_correlation_area_1 : %u:%u:%u\n', flat_corr_area1(1), flat_corr_area1(2) - flat_corr_area1(1), flat_corr_area1(end));
fprintf(fid, 'flat_field_correlation_area_2 : %u:%u:%u\n', flat_corr_area2(1), flat_corr_area2(2) - flat_corr_area2(1), flat_corr_area2(end));
if ~read_flatcor
    fprintf(fid, 'min_max_of_all_darks : %6g %6g\n', dark_min, dark_max);
    fprintf(fid, 'min_max_of_median_dark : %6g %6g\n', dark_med_min, dark_med_max);
    fprintf(fid, 'min_max_of_all_flats : %6g %6g\n', flat_min, flat_max);
    fprintf(fid, 'min_max_of_all_corrected_flats : %6g %6g\n', flat_min2, flat_max2);
    fprintf(fid, 'min_max_of_all_raws :  %6g %6g\n', raw_min, raw_max);
    fprintf(fid, 'min_max_of_all_corrected_raws :  %6g %6g\n', raw_min2, raw_max2);
    fprintf(fid, 'min_max_of_all_flat_corr_projs : %g %g \n', proj_min, proj_max);
end
% Phase retrieval
if do_phase_retrieval(1)
    fprintf(fid, 'do_phase_retrieval : %u\n', do_phase_retrieval);
    fprintf(fid, 'phase_retrieval_method : %s\n', phase_retrieval_method);
    fprintf(fid, 'phase_retrieval_regularisation_parameter : %f\n', phase_retrieval_reg_par);
    fprintf(fid, 'phase_retrieval_binary_filter_threshold : %f\n', phase_retrieval_bin_filt);
    fprintf(fid, 'phase_retrieval_cutoff_frequ : %f pi\n', phase_retrieval_cutoff_frequ / pi);
    fprintf(fid, 'phase_padding : %u\n', phase_padding);
end
% Volume
fprintf(fid, 'volume_shape : %u %u %u\n', vol_shape(1), vol_shape(2), vol_shape(3));
fprintf(fid, 'volume_size : %f %f %f %f %f %f\n', vol_size(1), vol_size(2), vol_size(3), vol_size(4), vol_size(5), vol_size(6));
% Rotation
fprintf(fid, 'excentric_rot_axis : %i\n', excentric_rot_axis);
fprintf(fid, 'crop_at_rot_axis : %u\n', crop_at_rot_axis);
fprintf(fid, 'stitch_projections : %u\n', stitch_projections);
fprintf(fid, 'stitch_method : %s\n', stitch_method );
fprintf(fid, 'rotation_angle_full_rad : %f * pi\n', rot_angle_full / pi);
fprintf(fid, 'rotation_angle_offset_rad : %f * pi\n', rot_angle_offset / pi);
fprintf(fid, 'rotation_axis_offset_calculated : %f\n', rot_axis_offset_calc);
fprintf(fid, 'rotation_axis_offset : %f\n', rot_axis_offset);
fprintf(fid, 'rotation_axis_position_calculated : %f\n', rot_axis_pos_calc);
fprintf(fid, 'rotation_axis_position : %f\n', rot_axis_pos);
fprintf(fid, 'rotation_axis_tilt_calculated : %f\n', rot_axis_tilt_calc);
fprintf(fid, 'rotation_axis_tilt : %f\n', rot_axis_tilt);
fprintf(fid, 'raw_image_binned_center : %f\n', raw_im_shape_binned1 / 2);
fprintf(fid, 'rotation_correlation_area_1 : %u:%u:%u\n', rot_corr_area1(1), rot_corr_area1(2) - rot_corr_area1(1), rot_corr_area1(end));
fprintf(fid, 'rotation_correlation_area_2 : %u:%u:%u\n', rot_corr_area2(1), rot_corr_area2(2) - rot_corr_area2(1), rot_corr_area2(end));
fprintf(fid, 'interactive_determination_of_rot_axis : %u\n', interactive_determination_of_rot_axis);
% FBP
fprintf(fid, 'ring_filter : %u\n', ring_filter);
fprintf(fid, 'ring_filter_method : %s\n', ring_filter_method);
fprintf(fid, 'ring_filter_median_width : %s\n', sprintf( '%u ', ring_filter_median_width) );
fprintf(fid, 'fbp_filter_type : %s\n', fbp_filter_type);
fprintf(fid, 'fpb_filter_freq_cutoff : %f\n', fpb_filter_freq_cutoff);
fprintf(fid, 'fbp_filter_padding : %u\n', fbp_filter_padding);
fprintf(fid, 'fbp_filter_padding_method : %s\n', fbp_filter_padding_method);
fprintf(fid, 'butterworth_filter : %u\n', butterworth_filter);
fprintf(fid, 'butterworth_order : %u\n', butterworth_order);
fprintf(fid, 'butterworth_cutoff_frequency : %f\n', butterworth_cutoff_frequ);
fprintf(fid, 'astra_pixel_size : %f\n', astra_pixel_size);
fprintf(fid, 'take_negative_logarithm : %u\n', take_neg_log);
fprintf(fid, 'gpu_name : %s\n', gpu.Name);
fprintf(fid, 'min_max_of_all_slices : %g %g\n', vol_min, vol_max);
fprintf(fid, 'write_float : %u\n', write_float);
fprintf(fid, 'write_16bit : %g\n', write_16bit);
fprintf(fid, 'write_8bit : %u\n', write_8bit);
fprintf(fid, 'write_8bit_binned : %u\n', write_8bit_binned);
fprintf(fid, 'compression_method : %s\n', compression);
fprintf(fid, 'compression_limits : %f %f\n', tlow, thigh);
fprintf(fid, 'reco_bin : %u\n', reco_bin);
fprintf(fid, 'full_reconstruction_time : %.1f s\n', toc);
fprintf(fid, 'date_of_reconstruction : %s', datetime);
fprintf(fid, 'rotation_axis_offset : %f\n', rot_axis_offset);
%fprintf(fid, ' : \n', );
fclose(fid);
PrintVerbose(verbose, '\n log file : %s', logfile_name)

PrintVerbose(verbose && interactive_determination_of_rot_axis, '\nTime elapsed in interactive mode: %g s (%.2f min)', tint, tint / 60 );
PrintVerbose(verbose, '\nTime elapsed for computation: %g s (%.2f min)', toc - tint, (toc - tint) / 60 );
PrintVerbose(verbose, '\nFINISHED: %s\n\n', scan_name)
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dbclear if error
