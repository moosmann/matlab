function p07_ivw_20210701_11012199_reco( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Created on 12-Jul-2021 by moosmanj

if nargin < 1
    SUBSETS = [];
end
if nargin < 2
    RUN_RECO = 0;
end
if nargin < 3
    PRINT_PARAMETERS = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameter here. Otherwise parameters are taken from the current
% version of 'p05_reco' if not set in the section below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%
par.scan_path = pwd; % string/pwd. pwd: change to directory of the scan to be reconstructed, string: absolute scan path
par.ref_path = {}; % cell of strings. Additonal data sets to be included for the correlation of projections and reference images
par.read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
par.read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
par.read_flatcor_trafo = @(im) im; %fliplr( im ); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
par.read_sino_folder = ''; % subfolder to scan path
par.read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.energy = []; % eV! if empty: read from log file (log file values can be ambiguous or even missing sometimes)
par.sample_detector_distance = []; % in m. if empty: read from log file
par.eff_pixel_size = []; %1.07e-6; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
par.pixel_scaling = []; % to account for mismatch of eff_pixel_size with, ONLY APPLIED BEFORE TOMOGRAPHIC RECONSTRUCTION, HAS TO BE CHANGED!
par.read_image_log = 0; % bool, default: 0. Read metadata from image log instead hdf5, if image log exists
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.raw_bin = 2; % projection binning factor: integer
par.raw_roi = []; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
par.im_trafo = '';% 'rot90(im,1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
% STITCHING/CROPPING only for scans without lateral movment. Legacy support
par.crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
par.stitch_projections = 0; % for 2 pi cans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
par.stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
par.proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
par.ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
par.wo_crop = 0; % Do not crop images to account for random lateral shift
par.virt_s_pos = 0; % Correct sample position in reconsructed volume if virtual sample position motors are used
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.02 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.02 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [5 5]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
par.ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'median';'none';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
image_correlation.force_calc = 0; % bool. force compuation of correlation even though a (previously computed) corrlation matrix exists
image_correlation.num_flats = 19; % number of best maching flat fields used for correction
image_correlation.area_width = [1 100];%[-100 1];% correlation area: index vector or relative/absolute position of [first pix, last pix], negative indexing is supported
image_correlation.area_height = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.filter = 1; % bool, filter ROI before correlation
image_correlation.filter_type = 'median'; % string. correlation ROI filter type, currently only 'median' is implemnted
image_correlation.filter_parameter = {[3 3], 'symmetric'}; % cell. filter paramaters to be parsed with {:}
% 'median' : using medfilt2, parameters: {[M N]-neighboorhood, 'padding'}
% 'wiener' : using wiender2, parameters: {[M N]-neighboorhood}
ring_filter.apply = 1; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 1; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft_dec_levels = 1:6; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft_wname = 'db7';'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft_sigma = 3; %  suppression factor for 'wavelet-fft'
ring_filter.jm_median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
par.strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
par.norm_sino = 0; % not recommended, can introduce severe artifacts, but sometimes improves quality
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'tieNLO_Schwinger';'dpc';'tie';'qp';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
% Interactive phase retrieval not supported for method 'tieNLO_Schwinger'
phase_retrieval.reg_par = 1.1; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
phase_retrieval.tieNLO_Schwinger.sn = 10; % Schwinger regularization: points of support
phase_retrieval.tieNLO_Schwinger.smax = 10; % Schwinger regularization: maximumg support range
phase_retrieval.dpc_steps = 5;
phase_retrieval.dpc_bin = 4;
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, use to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D';'slice';  % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = [];%[-1.5 1.5 -1.5 1.5 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_full_range = [];% 2 * pi ;[]; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle_offset = pi; % global rotation of reconstructed volume
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR posITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.flip_scan_position = 0; % for debugging
tomo.rot_axis_tilt_camera = 0; % in rad. camera tilt w.r.t rotation axis.
tomo.rot_axis_tilt_lamino = 0; % in rad. lamino tilt w.r.t beam.
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
tomo.take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
tomo.algorithm =  'fbp';'cgls';'sirt';'sart';'em';'fbp-astra'; % SART/EM only work for 3D reco mode
tomo.iterations = 50; % for iterateive algorithms: 'sirt', 'cgls', 'sart', 'em'
tomo.MinConstraint = []; % sirt3D/sirt2d/sart2d only. If specified, all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions, for example.
tomo.MaxConstraint = []; % sirt3D/sirt2d/sart2d only. If specified, all values above MaxConstraint will be set to MaxConstraint.
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID, mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = '';% sprintf( '%s%04u, tomo.algorithm, tomo.iterations); '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder_flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder_phase_map = ''; % subfolder in 'phase_map'
write.subfolder_sino = ''; % subfolder in 'sino'
write.subfolder_reco = ''; % subfolder in 'reco'
write.flatcor = 0; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression_method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression_method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 0; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression_method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression_method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression_method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression_parameter = [0.02 0.02]; % compression-method specific parameters
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 1; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 1; % find laminography tilt instead camera tilt
interactive_mode.angles = 0; % reconstruct slices with different scalings of angles
interactive_mode.angle_scaling_default_search_range = []; % if empty: use a variaton of -/+5 * (angle increment / maximum angle)
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 1; % Interactive retrieval to determine regularization parameter
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.poolsize = 0.6;
par.poolsize_gpu_limit_factor = 0.7; % Relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
par.gpu_index = tomo.astra_gpu_index;
par.use_gpu_in_parfor = 1;
phase_retrieval.use_parpool = 0; % bool. Disable parpool when out-of-memory error occurs during phase retrieval.
par.window_state = 'minimized';'normal';'maximized'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p07/2021/data/11012199/raw/';

image_correlation.num_flats = 19;
image_correlation.area_width = [301 400];
image_correlation.area_height = [0.2 0.8];
tomo.interpolate_missing_angles = 1;
interactive_mode.rot_axis_pos = 0;

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
tomo.rot_axis_offset = 1.5 * 2 / par.raw_bin;
par.scan_path = [raw_path 'ivw0015_referenzblau02c']; ADD

tomo.vol_size = [-0.5 0 -0.25 0.25 -0.5 0.5];
tomo.rot_axis_offset = -1.25 * 2 / par.raw_bin;
par.scan_path = [raw_path 'ivw0015_referenzblau_1_000']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_1_001']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_1_002']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_1_003']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_1_004']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_1_005']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_1_006']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_1_007']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_1_008']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
interactive_mode.rot_axis_pos = 0;
par.scan_path = [raw_path 'ivw0015_referenzblau_2_000']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_001']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_002']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_003']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_004']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_005']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_006']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_007']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_008']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_009']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_010']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_011']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_012']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_013']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_014']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_015']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_016']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2_017']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0015_referenzblau_2b_000']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2b_001']; ADD
par.scan_path = [raw_path 'ivw0015_referenzblau_2b_002']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
%% check rot axis
par.scan_path = [raw_path 'ivw0015_referenzblau_2c_000']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0016_referenzblau_3_000']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_001']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_002']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_003']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_004']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_005']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_006']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_007']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_008']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_009']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_010']; ADD
par.scan_path = [raw_path 'ivw0016_referenzblau_3_011']; ADD

tomo.vol_size = [-0.4 0.1 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_000']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_001']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_002']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_003']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_004']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_005']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_006']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_007']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_008']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_009']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_010']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_011']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_012']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_013']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_014']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_015']; ADD
par.scan_path = [raw_path 'ivw0017_Struktur1_gruen_1_016']; ADD

tomo.vol_size = [-0.4 0.1 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0018_Struktur1_gruen_1b_000']; ADD
par.scan_path = [raw_path 'ivw0018_Struktur1_gruen_1b_001']; ADD
par.scan_path = [raw_path 'ivw0018_Strukturierung1gruen_1b']; ADD

tomo.vol_size = [0 0.25 0 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0019_Strukturierung1gruen_2_nullscan']; ADD

tomo.vol_size = [-0.2 0.2 -0.2 0.2 -0.5 0.5];
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_000']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_001']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_002']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_003']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_004']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_005']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_006']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_007']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_008']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_009']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_010']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_011']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_012']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_013']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_014']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_015']; ADD
par.scan_path = [raw_path 'ivw0020_Struktur1_gruen_2_016']; ADD

tomo.vol_size = [-0.2 0.2 -0.2 0.2 -0.5 0.5];
par.scan_path = [raw_path 'ivw0021_Struktur1_gruen_2b_000']; ADD
par.scan_path = [raw_path 'ivw0021_Struktur1_gruen_2b_001']; ADD
par.scan_path = [raw_path 'ivw0021_Struktur1_gruen_2b_002']; ADD
par.scan_path = [raw_path 'ivw0021_Struktur1_gruen_2b_003']; ADD
par.scan_path = [raw_path 'ivw0021_Struktur1_gruen_2b_004']; ADD

tomo.vol_size = [-0.2 0.2 -0.2 0.2 -0.5 0.5];
par.scan_path = [raw_path 'ivw0022_Struktur1_gruen_2c_000']; ADD
par.scan_path = [raw_path 'ivw0022_Struktur1_gruen_2c_001']; ADD
par.scan_path = [raw_path 'ivw0022_Struktur1_gruen_2c_002']; ADD

tomo.vol_size = [-0.2 0.2 -0.3 0.1 -0.5 0.5];
par.scan_path = [raw_path 'ivw0023_Struktur1_gruen_3_000']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_000']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_001']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_002']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_003']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_004']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_005']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_006']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_007']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_008']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_009']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_010']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_011']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_012']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_013']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_014']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_015']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_016']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_017']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_018']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_019']; ADD
par.scan_path = [raw_path 'ivw0025_Struktur1_gruen_3_020']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_000']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_001']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_002']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_003']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_004']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_005']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_006']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_007']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_008']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_009']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_010']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_011']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_012']; ADD
par.scan_path = [raw_path 'ivw0026_Struktur2_pink_1_013']; ADD

tomo.vol_size = [-0.5 0 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_1b_000']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_1b_001']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_1b_002']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_1b_003']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_1b_004']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_1b_005']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_1b_006']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_000']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_001']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_002']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_003']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_004']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_005']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_006']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_007']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_008']; ADD
par.scan_path = [raw_path 'ivw0027_Struktur2_pink_2_009']; ADD

tomo.vol_size = [-0.3 0.2 -0.3  0.2 -0.5 0.5];
par.scan_path = [raw_path 'ivw0028_Struktur2_pink_1c_000']; ADD
par.scan_path = [raw_path 'ivw0028_Struktur2_pink_2c_000']; ADD

tomo.vol_size = [-0.2 0.2 -0.2  0.2 -0.5 0.5];
par.scan_path = [raw_path 'ivw0029_Struktur2_pink_2c_aperture05_000']; ADD

tomo.vol_size = [-0.2 0.2 -0.2  0.2 -0.5 0.5];
par.scan_path = [raw_path 'ivw0030_Struktur2_pink_2c_aperture07_000']; ADD

tomo.vol_size = [-0.2 0.2 -0.2  0.2 -0.5 0.5];
par.scan_path = [raw_path 'ivw0031_Struktur2_pink_2c_8001proj_000']; ADD

tomo.vol_size = [-0.5 0.0 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0032_Referenz_blau_4_000']; ADD

tomo.vol_size = [-0.5 0.0 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0033_Referenz_blau_5_000']; ADD
par.scan_path = [raw_path 'ivw0033_Referenz_blau_5_001']; ADD
par.scan_path = [raw_path 'ivw0033_Referenz_blau_5_002']; ADD
par.scan_path = [raw_path 'ivw0033_Referenz_blau_5_003']; ADD
par.scan_path = [raw_path 'ivw0033_Referenz_blau_5_004']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0034_Referenz_blau_6_35N']; ADD

tomo.vol_size = [-0.25 0.25 -0.25 0.25 -0.5 0.5];
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_000']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_001']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_002']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_003']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_004']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_005']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_006']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_007']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_008']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_009']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_010']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_011']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_012']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_013']; ADD
par.scan_path = [raw_path 'ivw0035_Struktur2_pink_3_014']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
