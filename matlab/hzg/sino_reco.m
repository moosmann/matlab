function sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range)
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

% if ~exist( 'external_parameter' ,'var')
%     clearvars
% end
close all hidden % close all open windows
dbstop if error

if nargin < 1
    scan_path = ... % pwd;
    '/asap3/petra3/gpfs/p07/2023/data/11017607/processed/bmc010_brainD_wholeEtOH_zM72p4_10rings';
    '/asap3/petra3/gpfs/p07/2023/data/11017607/processed/bmc003_brainB_slice4_paraffin_5p6mm_8rings';
    '/asap3/petra3/gpfs/p07/2023/data/11017598/processed/tig002_81217_mousebrain_202309_DESY';
    '/asap3/petra3/gpfs/p07/2023/data/11019133/processed/bmc009_81212_mousebrain_202311_DESY';
    '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw009_cet452a_FF99_Bp';
    '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw007_cet488b_F18096_Bp';
    '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw013_cet518a_MB7_Mb';
    '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw011_cet495b_M548_20_Ha';
    '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw012_cet548a_OO01_Oo';
    '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw006_cet547a_sd01_pp';
    '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw010_cet433b_UT1611_Pp';
    '/asap3/petra3/gpfs/p07/2023/data/11017607/processed/mosaic/bmc003_B4_paraffin/tests_hs01';
    '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw006_cet547a_sd01_pp';
    '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc015_B4_bobble_wt0p25_75px_b25_';
    '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2';
    '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc015_B4_bobble_wt0p25_75px_b25_';
    '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc003_B2_8rings';
    '/asap3/petra3/gpfs/p07/2022/data/11015573/processed/hereon_bmc_brain_small_roi_scan_44keV';
    '/asap3/petra3/gpfs/p07/2020/data/11010206/processed/mgbone05_19118_ti_8m';
    '/asap3/petra3/gpfs/p05/2020/data/11008823/processed/nova004_pyrochroa_coccinea_a';
    '/asap3/petra3/gpfs/p05/2020/data/11008823/processed/nova003_pentatomidae';
    '/asap3/petra3/gpfs/p07/2019/data/11006991/processed/hzg_ind_01_cork_a/';
    '/asap3/petra3/gpfs/p05/2020/data/11010107/processed/bmc05_v63l';
    '/asap3/petra3/gpfs/p05/2020/data/11010107/processed/bmc07_v67r';
    '/asap3/petra3/gpfs/p05/2019/data/11007580/processed/smf_09_be_3033';
    '/asap3/petra3/gpfs/p07/2019/data/11007454/processed/bmc06_tooth1';
    '/asap3/petra3/gpfs/p07/2022/data/11015567/processed/mbs020_sample1_mgti_a';
end
if nargin < 2
    rot_angle_full_range = [];
end
if nargin < 3
    rot_axis_search_range = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scan_path = '';pwd;
raw_bin = 2; % projection binning factor: integer
proj_range = []; % projection range
%read_sino_folder = sprintf( 'trans%02u', raw_bin);% string. default: '', picks first trans folder found
read_sino_folder = sprintf( 'trans%02u_360', raw_bin);% string. default: '', picks first trans folder found
%read_sino_folder = sprintf( 'test04', raw_bin);% string. default: '', picks first trans folder found
%read_sino = 1; % read preprocessed sinograms. CHECK if negative log has to be taken!
read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
sino_range = [];
energy = []; % in eV!
sample_detector_distance = []; % in m
eff_pixel_size = []; % in m, read from reconlog.txt if []
%%% RING FILTER
ring_filter.apply = 1; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft_dec_levels = 1:6; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft_wname = 'db7';'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft_sigma = 3; %  suppression factor for 'wavelet-fft'
ring_filter.jm_median_width = [3 7 11]; % multiple widths are applied consecutively, eg [3 11 21 31 39];
%%% PIXEL FILTER
filter_sino = 1; % Filter hot pixesl
pixel_filter_threshold = [0.001 0.00005]; % threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [3 3];
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
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
tomo.reco_mode = 'slice';%'3D'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.slab_wise = 1;
tomo.slices_per_slab = [];
tomo.vol_size = [];%[-.1 0.1 -0.1 0.1 -0.5 0.5];%[-1 1 -1 1 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_offset = 0; % global rotation of reconstructed volume
tomo.rot_angle_full_range = rot_angle_full_range; % in rad, read from reconlog if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_axis_offset = 0;%[]; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.rot_axis_tilt_camera = 0; % in rad. camera tilt w.r.t rotation axis.
%tomo.rot_axis_tilt_camera = 0.0028;
tomo.rot_axis_tilt_lamino = 0; %
tomo.rot_axis_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
tomo.rot_axis_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
tomo.fbp_filter_type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
tomo.fbp_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
tomo.fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
tomo.fbp_filter_padding_method = 'symmetric';
tomo.butterworth_filter = 0; % use butterworth filter in addition to FBP filter
tomo.butterworth_filter_order = 1;
tomo.butterworth_filter_frequ_cutoff = 0.9;
tomo.astra_pixel_size = 1; % detector pixel size for reconstruction: if different from one 'tomo.vol_size' must to be ajusted, too!
tomo.take_neg_log = 1; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
tomo.algorithm = 'fbp';'sirt'; 'cgls';'sart';'em';'fbp-astra'; % SART/EM only work for 3D reco mode
tomo.iterations = 40; % for 'sirt' or 'cgls'.
tomo.sirt_MinConstraint = []; % If specified, all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions, for example.
tomo.sirt_MaxConstraint = []; % If specified, all values above MaxConstraint will be set to MaxConstraint.
tomo.rot_axis_search_auto = 0; % find extrema of metric within search range
tomo.rot_axis_search_range = rot_axis_search_range; % search reach for automatic determination of the rotation axis offset, overwrite interactive result if not empty
tomo.rot_axis_search_metric = 'neg'; % string: 'neg','entropy','iso-grad','laplacian','entropy-ML','abs'. Metric to find rotation axis offset
tomo.rot_axis_search_extrema = 'max'; % string: 'min'/'max'. chose min or maximum position
tomo.rot_axis_search_fit = 1; %
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; % absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID, mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder_flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder_phase_map = ''; % subfolder in 'phase_map'
write.subfolder_sino = ''; % subfolder in 'sino'
write.subfolder_reco = '';%read_sino_folder;''; % subfolder in 'reco'
write.flatcor = 0; % save preprocessed flat corrected projections
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression_method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression_method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 1; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression_method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression_method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression_method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression_parameter = [0.02 0.02]; % compression-method specific parameter
write.outputformat = 'tif';'hdf_volume';'hdf_slice'; % string
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 0; % show images and plots during reconstruction
interactive_mode.show_stack_imagej = 1; % use imagej instead of MATLAB to scroll through images during interactive mode
interactive_mode.show_stack_imagej_use_virtual = 0; % use virtual stack for faster loading, but slower scrolling
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
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
par.use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.poolsize = 0.5; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.5; % Relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
par.use_gpu_in_parfor = 1;
par.gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], use all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
logpar = [];
if interactive_mode.rot_axis_tilt && strcmpi( tomo.reco_mode, 'slice' )
    error( 'Slicewise reconstruction and reconstruction with tilted rotation axis are not compatible!' )
end
% infer (relative) vol_shape from vol_size if vol_shape is empty, but vol_size is given
if ~isempty( tomo.vol_size ) && isempty( tomo.vol_shape )
    tomo.vol_shape = tomo.vol_size(2:2:end) - tomo.vol_size(1:2:end);
end
if ~isempty( tomo.rot_axis_offset ) && ~isempty( tomo.rot_axis_position )
    error('tomo.rot_axis_offset (%f) and tomo.rot_axis_position (%f) cannot be used simultaneously. One must be empty.', tomo.rot_axis_offset, tomo.rot_axis_position)
end

% Default assignment if non-existing or empty!
assign_default( 'tomo.rot_axis_offset', [] )
assign_default( 'pixel_filter_radius', [3 3] )
assign_default( 'image_correlation.force_calc', 0 );
assign_default( 'write.path', '' )
assign_default( 'write.parfolder', '' )
assign_default( 'write.subfolder_reco', '' )
assign_default( 'write.subfolder_flatcor', '' )
assign_default( 'write.subfolder_phase_map', '' )
assign_default( 'write.subfolder_sino', '' )
assign_default( 'write.sino_shift_cropped', 0 )
assign_default( 'write.deleteFiles', 0)
assign_default( 'write.beamtimeID', '' )
assign_default( 'tomo.reco_mode', '3D' )
assign_default( 'tomo.rot_axis_tilt', 0 )
assign_default( 'write.scan_name_appendix', '' )
assign_default( 'interactive_mode.rot_axis_pos_default_search_range', -4:0.5:4 ) % binned pixel
assign_default( 'interactive_mode.rot_axis_tilt_default_search_range', -0.005:0.001:0.005 ) % radian
assign_default( 'interactive_mode.phase_retrieval_default_search_range', [] )
assign_default( 'interactive_mode.angles', 0 );
assign_default( 'interactive_mode.angle_scaling_default_search_range', [] );
assign_default( 'tomo.rot_axis_corr_area2', [0.1 0.9] );
assign_default( 'par.window_state', 'normal' );
assign_default( 'tomo.rot_axis_search_auto', 0);
assign_default( 'tomo.rot_axis_search_metric', '');
assign_default( 'tomo.rot_axis_search_verbose', 1);
assign_default( 'tomo.rot_axis_search_extrema', 'max' );
assign_default( 'tomo.rot_axis_search_auto', 1);
assign_default( 'par.gpu_index', [])
assign_default( 'write.outputformat', 'tif')

% Define variables from struct fields for convenience
raw_bin = single( raw_bin );
par.raw_bin = raw_bin;
eff_pixel_size_binned = raw_bin * eff_pixel_size;
outputformat = write.outputformat;
par.verbose = 0;

astra_clear % if reco was aborted, ASTRA memory is not cleared

% Utility functions
%imsc1 = @(im) imsc( rot90( im ) );
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
write.scan_name = scan_name;
[beamtime_path, ~] = fileparts(raw_path);
[~, beamtime_id] = fileparts(beamtime_path);

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

% interactive path
p = [write.path, filesep, 'interactive', filesep];
p = regexprep( p, 'processed', 'scratch_cc');
write.interactive_path = p;

% Save raw path to file for shell short cut
filename = [userpath, filesep, 'path_to_raw'];
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

% Start parallel CPU pool %%%
%if ~strcmp(tomo.reco_mode,'slice') || interactive_mode.rot_axis_pos == 1
[~, par.poolsize] = OpenParpool( par.poolsize , par.use_cluster, [beamtime_path filesep 'scratch_cc']);
% GPU info
if isempty( par.gpu_index )
    par.gpu_index = 1:gpuDeviceCount;
end
tomo.gpu_index = par.gpu_index;
% GPU info quick
% fprintf( '\n GPU : [index, total memory/GiB] =' )
% num_gpu = numel( par.gpu_index );
% for mm = 1:num_gpu
%     nn = par.gpu_index(mm);
%     gpu = parallel.gpu.GPUDevice.getDevice( nn );
%     mem_total = gpu.TotalMemory/1024^3;
%     fprintf( ' [%u %.3g]', nn, mem_total )
% end
% GPU info detailed, needed for parpool optimization
mem_avail_gpu = zeros([1, numel(par.gpu_index)]);
mem_total_gpu = mem_avail_gpu;
gpu_index = par.gpu_index;
num_gpu = numel( par.gpu_index );

parfor mm = 1:num_gpu
    nn = gpu_index(mm);        
    gpu = gpuDevice(nn);
    gpu.reset;
    mem_avail_gpu(mm) = gpu.AvailableMemory;
    mem_total_gpu(mm) = gpu.TotalMemory;
    ma = mem_avail_gpu(mm)/1024^3;
    mt = mem_total_gpu(mm)/1024^3;
    r = 100 * ma / mt;
    fprintf( '\n GPU %u %u: memory: total: %.3g GiB, available: %.3g GiB (%.2f%%)', mm,nn, mt, ma, r)
end
par.mem_avail_gpu = mem_avail_gpu;
par.mem_total_gpu = mem_total_gpu;
%end

% Save scan path to file
filename = [userpath, filesep, 'path_to_scan'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', scan_path );
fclose( fid );

% Reco path
if isempty( write.subfolder_reco )
    reco_path = [write.path, filesep, 'reco', filesep];
else
    reco_path = [write.path, filesep, 'reco', filesep, write.subfolder_reco, filesep];
end
write.reco_path = reco_path;

% Figure path
write.fig_path = [write.path, filesep, 'figures', filesep];
fig_path = write.fig_path;
CheckAndMakePath( fig_path )

% read_image_par = @(filename) read_image( filename, par );
filt_pix_par.threshold_hot = pixel_filter_threshold(1);
filt_pix_par.threshold_dark = pixel_filter_threshold(2);
filt_pix_par.medfilt_neighboorhood = pixel_filter_radius;
filt_pix_par.filter_dead_pixel = 1;
filt_pix_par.filter_Inf = 0;
filt_pix_par.filter_NaN = 1;
filt_pix_par.verbose = 0;
filt_pix_par.use_gpu = par.use_gpu_in_parfor;

t = toc;
%% Read sinogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sino path
if isempty( read_sino_folder )
    d = dir([scan_path filesep 'trans*']);
    if numel( d ) < 1
        error( '''trans'' folder not found. Check ''scan_path'' and ''read_sino_folder''.' )
    else
        fprintf( '\n %u trans folder(s) found. Set ''read_sino_folder'' to ''%s''', numel(d ), d(1).name )
        sino_path = [scan_path d(1).name filesep];
    end
else
     sino_path = [scan_path read_sino_folder filesep];
end

% Sinogram file names
data_struct = dir( [sino_path filesep '*.tif'] );
if isempty( data_struct )
    error('Sinograms not found.')
else
    sino_names = {data_struct.name};
end
if ~isempty(sino_range)
    if numel( sino_range ) == 2
        sino_range = sino_range(1):sino_range(2);
    end
    sino_names = sino_names(sino_range);
end
im_shape_binned2 = numel( sino_names );
sino_names_mat = NameCellToMat( sino_names );

% Parameters
filename = sprintf('%s%s', sino_path, sino_names_mat( round( im_shape_binned2 / 2 ), :));
par.tifftrafo = 0;
[sino, ~] = read_image( filename, par );
sino = read_sino_trafo( sino );
im_shape_cropbin1 = size( sino, 2 );
%num_proj_read = size( sino, 1 );
if ~isempty( proj_range )
    sino = sino(proj_range, :);
end
num_proj_used = size( sino, 1 );

% Read parameters from reconlog.txt
fn = [scan_path filesep 'reconlog.txt'];
if exist(fn, 'File')
    fprintf('\nReading reconlog.txt:')
    fid = fopen(fn);
    c = textscan(fid, '%s : %s');
    fclose(fid);
    
    % angle increment
    b = contains(c{1},'scan_angularstep');
    if sum(b)
        angularstep = c{2}(b);
        angularstep = pi / 180 *str2double(angularstep{1});
        fprintf('\n angle increment : %f rad = %f deg', angularstep, angularstep * 180 / pi)
    end
    
    % number of projections
    b = contains(c{1},'scan_n_angle');
    if sum(b)
        n_angle = c{2}(b);
        n_angle = str2double(n_angle{1});
        fprintf('\n number of angles : %u', n_angle)
        angles = (0:n_angle - 1) * angularstep;
        ar = n_angle * angularstep;
        if isempty(tomo.rot_angle_full_range)
            tomo.rot_angle_full_range = ar;
        end
    end
    
    % pixel size
    b = contains(c{1},'scan_pixelsize');
    if sum(b)
        pixelsize = c{2}(b);
        pixelsize = str2double(pixelsize{1});
        fprintf('\n pixelsize : %f', pixelsize)
        c{1}(b) = [];
        c{2}(b) = [];
    end
    
    b = contains(c{1},'pixelsize');
    if sum(b)
        pixelsize = str2double(c{2}{b});
        fprintf('\n pixelsize : %f micron', pixelsize * 1000)
        eff_pixel_size = raw_bin * pixelsize / 1000;
    end
    
    b = contains(c{1},'reco_mode');
    if sum(b)
        reco_mode = str2double(c{2}{b});
        fprintf('\n reco mode : %f', reco_mode)
        ar =  reco_mode * pi / 180;
    end
    
    if isempty(tomo.rot_angle_full_range)
        tomo.rot_angle_full_range = ar;
    end
    
    if isempty(eff_pixel_size)
        eff_pixel_size = raw_bin * pixelsize;
    end
end
ar = tomo.rot_angle_full_range;
fprintf('\n angular range : %f rad = %f * pi rad = %f deg', ar, ar / pi, ar * 180 / pi)
% Angles
if  ~exist( 'angles', 'var' )
    if isempty( tomo.rot_angle_full_range )
        cprintf( 'Red', '\nEnter full angle of rotation (including one additional increment) or vector of angles, in radians: ' );
        tomo.rot_angle_full_range = input( '' );
    end
    if isscalar( tomo.rot_angle_full_range )
        angles = tomo.rot_angle_full_range * (0:num_proj_used - 1) / num_proj_used;
    else
        angles = tomo.rot_angle_full_range;
    end
    if length( angles ) ~= num_proj_used
        error( 'Number of angles (%u) entered not consistent with sinogram (%u) read.', numel( angles), num_proj_used )
    end
end
if tomo.rot_angle_full_range > 90
    error('rotation angle full range, %f, is most likely given in degrees instead of rad', tomo.rot_angle_full_range)
end

nn = round(im_shape_binned2 / 2 );
filename = sprintf('%s%s', sino_path, sino_names_mat(nn,:));
fprintf('\n sino for interactive mode:\n %s',filename)

sino = read_image( filename );
sino = read_sino_trafo( sino );
if ~isempty( proj_range )
    sino = sino(proj_range, :);
end
%[s1, s2] = size( sino );
%sino = FilterPixel_par(sino);
proj = reshape( sino, [im_shape_cropbin1, 1, num_proj_used] );

%% TOMOGRAPHY: interactive mode to find rotation axis offset and tilt %%%%%
[tomo, angles, tint] = interactive_mode_rot_axis( par, logpar, phase_retrieval, tomo, write, interactive_mode, proj, angles);
tomo.angles = angles;

%% Automatic rot axis determination
tomo = find_rot_axis_offset_auto(tomo, proj, par, write, interactive_mode);

%% Tomographic reconstruction
fprintf( '\nTomographic reconstruction:')
fprintf( '\n method : %s', tomo.algorithm )
fprintf( '\n angle scaling : %g', tomo.angle_scaling )
fprintf( '\n angles [first last]/pi : [%g %g]', angles( [1 end] ) / pi )
if isscalar( tomo.rot_angle_full_range)
    ar = tomo.rot_angle_full_range;
else
    ar = max(tomo.rot_angle_full_range) - min(tomo.rot_angle_full_range);
end
fprintf('\n rotation angle range: %f * pi = %f deg', ar / pi, ar / pi * 180)
fprintf( '\n volume shape : [%g, %g, %g]', tomo.vol_shape )
%tomo.angles = tomo.rot_angle_offset + angles;

if write.reco
    %reco_bin = write.reco_binning_factor; % alias for readablity
    CheckAndMakePath( reco_path, 0 )
    % Save reco path to file
    filename = [userpath, filesep, 'path_to_reco'];
    fid = fopen( filename , 'w' );
    fprintf( fid, '%s', reco_path );
    fclose( fid );
end

filt = iradonDesignFilter(tomo.fbp_filter_type, (1 + tomo.fbp_filter_padding) * size( proj, 1), tomo.fbp_filter_freq_cutoff);
padding = tomo.fbp_filter_padding;
padding_method = tomo.fbp_filter_padding_method;
save_path = sprintf( '%sfloat_rawBin%u/%s', reco_path, raw_bin);
CheckAndMakePath( save_path );
fprintf('\n save path: %s',save_path)
take_neg_log = tomo.take_neg_log;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.verbose = 0;
FilterPixel_par = @(sino,gpu_index) FilterPixel(sino,filt_pix_par,gpu_index);
pp_filter_ring_artefacts_par = @(sino) pp_filter_ring_artefacts( ring_filter, sino, angles, par);
ring_filter_apply = ring_filter.apply;
if tomo.slab_wise
    fprintf( '\nStart slab-wise reconstruction')
    if isempty( tomo.slices_per_slab )
        num_slices_per_slab = 1 * par.poolsize;
        tomo.slices_per_slab = num_slices_per_slab;
    else
        num_slices_per_slab = tomo.slices_per_slab;
    end
    num_slabs = ceil( im_shape_binned2 / num_slices_per_slab );
    proj = zeros( [im_shape_cropbin1, num_slices_per_slab, num_proj_used], 'single' );
    fprintf( ' \n number of slices in total : %u', im_shape_binned2 )
    fprintf( ' \n number of slices per slab : %u', num_slices_per_slab )
    fprintf( ' \n number of slabs : %u', num_slabs )
    % Loop over slab
    % Hack to continue loop
    if exist('ll', 'var')
        start_slab = ll;
    else
        start_slab = 1;
    end
    %% Slab loop
    for ll = start_slab:num_slabs
        
        % Slab indices
        s0 = 1 + (ll - 1) * num_slices_per_slab;
        s1 = min( [num_slices_per_slab + (ll - 1) * num_slices_per_slab, im_shape_binned2] );
        slab_sino_names_mat = sino_names_mat(s0:s1,:);
        num_slices = s1 - s0 + 1;
        
        fprintf( '\nSlab no. %4u of %u', ll, num_slabs )
        fprintf( ', slab range : [%6u %6u]', s0, s1 )
        fprintf( ', num slices : %4u', num_slices )
        
        % Read and filter
        parfor mm = 1:num_slices
            % Read
            sino_name = slab_sino_names_mat(mm, :);
            filename = sprintf('%s%s', sino_path, sino_name);
            sino = read_image( filename );
            sino = read_sino_trafo( sino );
            if filter_sino
                gpu_index = mod(mm, num_gpu) + 1;
                sino = FilterPixel_par(sino,gpu_index);
            end
            if ~isempty( proj_range )
                sino = sino(proj_range,:,:);
            end
            % Check for errors
            num_zeros = sum( sino(:) == 0);
            num_nan = sum(isnan(sino(:)));
            num_inf = sum( isinf( sino(:)));
            if (num_zeros + num_nan + num_inf) > 0
                fprintf( '\n [%u %u] %s: %u zeros, %u NaN, %u Inf', ll, mm, sino_name, num_zeros, num_inf, num_nan)
            else
                sino = NegLog( sino, take_neg_log );
                sino = reshape( sino, [im_shape_cropbin1, 1, num_proj_used] );                
                if ring_filter_apply
                    sino = pp_filter_ring_artefacts_par(sino);
                end
                % Filter
                sino = padarray( sino, padding * [im_shape_cropbin1 0 0], padding_method, 'post' );
                sino = real( ifft( bsxfun(@times, fft( sino, [], 1), filt), [], 1, 'symmetric') );
                sino = sino(1:im_shape_cropbin1,:,:);
                proj(:,mm,:) = sino;
            end
        end
        
        % Tomo
        % Delete empty slices
        proj(:,num_slices + 1:end,:) = [];
        
        %% ADJUST tomo struct: volume size and shape
        [tomo.vol_shape, tomo.vol_size] = volshape_volsize( proj, [], [] );
        %[tomo.vol_shape, tomo.vol_size] = volshape_volsize( proj, tomo.vol_shape, tomo.vol_size, tomo.rot_axis_offset, 0 );
        vol = astra_parallel3D( tomo, permute( proj, [1 3 2]) );
        
        % Save
        switch outputformat
            case 'tif'
                parfor mm = 1:num_slices
                    sino_name = slab_sino_names_mat(mm,:);
                    filename = sprintf( '%s%s', save_path, sino_name);
                    write32bitTIFfromSingle(filename,vol(:,:,mm))
                end
            case 'hdf_slice'
                parfor mm = 1:num_slices
                    sino_name = slab_sino_names_mat(mm,:);
                    filename = sprintf( '%s%s.h5', save_path, sino_name(1:end-3));
                    write_hdf(filename,vol(:,:,mm),'slice')
                end
            case 'hdf_volume'
                [s1,s2,~] = size(vol);
                filename = sprintf( '%s%s.h5', save_path, scan_name);
                if ll == 1
                    h5create(filename,'/volume',[s1 s2 im_shape_binned2],'Datatype','single')
                end
                parfor mm = 1:num_slices
                    start = [1 1 mm + (ll - 1)* num_slices];
                    count = [s1 s2 1];
                    h5write(filename,'/volume',vol(:,:,mm),start,count)
                end
        end
    end
    
else
    fprintf( '\nStart slice-wise reconstruction:')
    % Slice-wise reco
    % Read sinogram, reco, write
    num_slices = numel( sino_names );
    num_gpu = gpuDeviceCount;
    parfor (nn = 1:num_slices,num_gpu)
        filename = sprintf('%s%s', sino_path, sino_names_mat(nn, :));
        gpu_index = mod(nn,num_gpu) + 1;
        sino = read_image( filename );
        sino = read_sino_trafo( sino );
        if filter_sino
            sino = FilterPixel(sino,filt_pix_par,gpu_index);
        end
        if ~isempty( proj_range )
            sino = sino(proj_range, :);
        end
        if sum( sino(:) == 0) || sum( isnan( sino(:) ) ) || sum( isinf( sino(:) ) )
            fprintf( ' %u', nn )
        else
            sino = NegLog( sino, take_neg_log );
            sino = reshape( sino, [im_shape_cropbin1, 1, num_proj_used] );
            % Ring filter
            if ring_filter_apply
                sino = pp_filter_ring_artefacts_par(sino);
            end
            % FBP Filter
            sino = padarray( sino, padding * [im_shape_cropbin1 0 0], padding_method, 'post' );
            sino = real( ifft( bsxfun(@times, fft( sino, [], 1), filt), [], 1, 'symmetric') );
            sino = sino(1:im_shape_cropbin1,:,:);
            % Tomo
            vol = astra_parallel2D( tomo, permute( sino, [3 1 2]), gpu_index);
            % Save
            switch outputformat
                case 'tif'
                    filename = sprintf( '%s/%s', save_path, sino_names_mat(nn,:));
                    write32bitTIFfromSingle(filename,vol)
                case 'hdf_slice'
                    filename = sprintf( '%s/%s.h5', save_path, sino_names_mat(nn,1:end-3));
                    write_hdf(filename,vol(:,:,nn),'slice')
                case 'hdf_volume'
                    filename = sprintf( '%s%s.h5', save_path, scan_name)
                    [s1,s2] = size(vol);
                    if ~exist(filename,'file')
                        h5create(filename,'/volume',[s1 s2 num_slices],'Datatype','single')
                    end
                    start = [1 1 nn];
                    count = [s1 s2 1];
                    h5write(filename,'/volume',vol,start,count)
            end
        end
    end
end
fprintf( ' \n tomo reco done in %.1f s (%.2f min)', toc-t, (toc-t)/60 )

%% Write reco log file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CheckAndMakePath( reco_path )
if write.reco
    logfile_path = reco_path;
end
CheckAndMakePath( logfile_path )
if write.reco
    logfile_name = sprintf( '%sreco_rawBin%u.log', logfile_path, raw_bin);
    
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
    if tomo.run
        % Volume
        fprintf(fid, 'tomo.vol_shape : %u %u %u\n', tomo.vol_shape(1), tomo.vol_shape(2), tomo.vol_shape(3));
        fprintf(fid, 'tomo.vol_size : %f %f %f %f %f %f\n', tomo.vol_size(1), tomo.vol_size(2), tomo.vol_size(3), tomo.vol_size(4), tomo.vol_size(5), tomo.vol_size(6));
        % Rotation
        fprintf(fid, 'tomo.reco_mode : %s\n', tomo.reco_mode);
        fprintf(fid, 'tomo.rot_angle_full_range : %f * pi rad\n', tomo.rot_angle_full_range / pi);
        fprintf(fid, 'tomo.rot_angle_offset : %f * pi rad\n', tomo.rot_angle_offset / pi);
        fprintf(fid, 'tomo.rot_axis_position : %f\n', tomo.rot_axis_position);
        fprintf(fid, 'tomo.rot_axis_tilt : %f\n', tomo.rot_axis_tilt);
        fprintf(fid, 'raw_image_binned_center : %f\n', im_shape_cropbin1 / 2);
        fprintf(fid, 'interactive_mode.rot_axis_pos : %u\n', interactive_mode.rot_axis_pos);
        fprintf(fid, 'interactive_mode.phase_retrieval : %u\n', interactive_mode.phase_retrieval);
        % Tomo
        fprintf(fid, 'tomo.algorithm : %s\n', tomo.algorithm );
        switch tomo.algorithm
            case 'fbp' %'fbp-astra'}
                fprintf(fid, 'tomo.fbp_filter_type : %s\n', tomo.fbp_filter_type);
                fprintf(fid, 'tomo.fbp_filter_freq_cutoff : %f\n', tomo.fbp_filter_freq_cutoff);
                fprintf(fid, 'tomo.fbp_filter_padding : %u\n', tomo.fbp_filter_padding);
                fprintf(fid, 'tomo.fbp_filter_padding_method : %s\n', tomo.fbp_filter_padding_method);
            case 'sirt'
                fprintf(fid, 'tomo.iterations : %u\n',  tomo.iterations );
                fprintf(fid, 'tomo.sirt_MinConstraint : %f\n',  tomo.sirt_MinConstraint);
                fprintf(fid, 'tomo.sirt_MaxConstraint : %f\n',  tomo.sirt_MaxConstraint);
            case 'cgls';'sart';'em';
                fprintf(fid, 'tomo.iterations : %u\n',  tomo.iterations );
        end
        fprintf(fid, 'tomo.butterworth_filter : %u\n', tomo.butterworth_filter);
        fprintf(fid, 'tomo.butterworth_filter_order : %u\n', tomo.butterworth_filter_order);
        fprintf(fid, 'tomo.butterworth_filter_frequ_cutoff : %f\n', tomo.butterworth_filter_frequ_cutoff);
        fprintf(fid, 'tomo.astra_pixel_size : %f\n', tomo.astra_pixel_size);
        fprintf(fid, 'tomo.take_neg_log : %u\n', tomo.take_neg_log);
        fprintf(fid, 'gpu_name : %s\n', gpuDevice().Name);
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
                fprintf(fid, 'write.compression_method : %s\n', write.compression_method);
            end
            if exist( 'tlow', 'var' ) && exist( 'thigh', 'var' )
                fprintf(fid, 'compression_limits : %f %f\n', tlow, thigh);
            end
            fprintf(fid, 'reco_bin : %u\n', write.reco_binning_factor);
        end
        fprintf(fid, 'full_reconstruction_time : %.1f s\n', toc);
        fprintf(fid, 'date_of_reconstruction : %s\n', datetime);
        fprintf(fid, 'tomo.rot_axis_offset at %u x binning : %f\n', raw_bin, tomo.rot_axis_offset);
    end
    fclose(fid);
    % End of log file
    fprintf( '\n log file : \n%s', logfile_name)
    fprintf( '\n reco_path : \n%s', reco_path)
end
PrintVerbose( interactive_mode.rot_axis_pos, '\nTime elapsed in interactive rotation axis centering mode: %g s (%.2f min)', tint, tint / 60 );
fprintf( '\nTime elapsed for computation: %g s (%.2f min)', toc - tint, (toc - tint) / 60 );
fprintf( '\nFINISHED: %s at %s\n', scan_name, datetime )
if exist('fast_reco','var') && fast_reco.run
    cprintf( 'Red', '\nATTENTION: fast reco mode was turned on!\n' )
end
fprintf( '\n')
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dbclear if error
