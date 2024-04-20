function p07_reco_loop_mousebrain_kuo_11017598( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Template function to loop over data sets given in the 'PARAMETER / DATA
% SETS' section below. The 'DEFAULT PARAMETERS' section defines the default
% paramters. To add a data / parameter to the loop, define your
% reconstruction parameters followed by 'ADD'. Using 'ADD('r') restores the
% default parameters after a datat set is added. Default parameters are
% defined with 'SET_DEFAULT' (or 'ADD_DEFAULT' or 'ADD('d')' or
% 'ADD('default')'.
%
% Caution: if parameters are changed, they remain changed after a data set
% was added unless 'ADD('r')' ('ADD('restore')) is used which restores 
% all parameters as defined in the 'DEFAULT PARAMETER' section by
% 'SET_DEFAULT'.
%
% Caution: all parameters not defined within the file are taken from main
% reconstruction script 'p05_reco'.
%
% Caution: if default parameters are not defined by closing the DEFAULT
% PARAMETER section with a 'SET_DEFAULT' statement, then the first call of
% 'ADD' will define default parameters.
% 
% Workflow:
% - Copy this file.
% - Check, modify, or copy parameter from 'p05_reco' in  'DEFAULT
%   PARAMETERS' section.
% - Add parameter / data sets to loop over in the ''PARAMETER / DATA SETS'
%   section. 
% - Type 'F5' or call without arguments to list all added data sets.
% - Choose subset of added data sets to loop over by indices (see SUBSETS
%   arguments below).
% - Use RUN_RECO equals 0 with PRINT_PARAMETERS (see below) to check
%   parameters setting.
% - Set RUN_RECO equals 1 to start the loop.
% 
% ARGUMENTS
% SUBSETS : 1D array of integers. subset of added data sets to be looped over
% RUN_RECO : bool. default: 0. 0: loops over the subsets but does not start
% reconstructions, 1: start the reconstruction loop.
% PRINT_PARAMETERS : string or cell of strings, eg {'raw_roi',
% 'parfolder'}. Parameter setting to be printed at each loop iteration.
% Useful in combination with RUN_RECO = 0 to check parameter setting for
% the subset to loop over
%
% Created on 07-Mar-2024 by moosmanj

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

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.scan_path = pwd; % string/pwd. pwd: change to directory of the scan to be reconstructed, string: absolute scan path, last_folder_modified('folder')
par.ref_path = {}; % cell of strings. Additonal data sets to be included for the correlation of projections and reference images
par.read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
par.read_flatcor_path = ''; % absolute path containing flat-field corrected projections
par.read_flatcor_range = 1; % scalar or vector. range of flatcorrected projections to be read
par.read_flatcor_bin = 1; % Binning of flat-corrected projections
par.read_flatcor_trafo = @(im) im; %fliplr( im ); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
par.read_sino_folder = ''; % subfolder to scan path
par.read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.read_sino_range = 1;
par.sino_roi = []; % horizontal ROI when reading sinograms, vertical ROI not yet implemented
par.filter_sino = 1; % bool. Pixel filtering of sinogram using parameters below:
pixel_filter_sino.threshold_hot = 0;
pixel_filter_sino.threshold_dark = 0;
pixel_filter_sino.medfilt_neighboorhood = [3 3];
pixel_filter_sino.filter_dead_pixel = 1;
pixel_filter_sino.filter_Inf = 1;
pixel_filter_sino.filter_NaN = 1;
pixel_filter_sino.verbose = 0;
par.energy = []; % eV! if empty: read from log file (log file values can be ambiguous or even missing sometimes)
par.sample_detector_distance = []; % in m. if empty: read from log file
par.eff_pixel_size = []; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
par.pixel_scaling = []; % to account for mismatch of eff_pixel_size with, ONLY APPLIED BEFORE TOMOGRAPHIC RECONSTRUCTION, HAS TO BE CHANGED!
par.read_image_log = 0; % bool, default: 0. Read metadata from image log instead hdf5, if image log exists
par.read_filenames_from_disk = 0; % only for stepscans with tiff subfolders
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.raw_bin = 2; % projection binning factor: integer
par.raw_roi = []; % vertical and/or horizontal ROI; coordinate (1,1) = top left pixel; supports absolute, relative, negative, and mixed indexing.
% []: use full image;
% [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2); if raw_roi(2) < 0 reads until end - |raw_roi(2)|; relative indexing similar.
% [y0 y1 x0 x1]: vertical + horzontal ROI, each ROI as above
% if scalar = 0: auto roi, selects vertical ROI automatically. Use only for DCM. Not working for *.raw data where images are flipped and DMM data.
% if scalar < 1: Threshold is set as min(proj(:,:,[1 end])) + abs(raw_roi)*median(dark(:)). raw_roi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
par.im_trafo = '';% 'rot90(im,1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
par.filter_ref = @(x) (x);
par.filter_proj = @(x) (x);
% STITCHING/CROPPING only for scans without lateral movment. Legacy support
par.crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
par.stitch_projections = 0; % for 2 pi cans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
par.stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
par.proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
par.ref_range = [];% range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
par.crop_proj = 0; % Crop images to account for random lateral shift
par.virt_s_pos = 0; % Correct sample position in reconsructed volume if virtual sample position motors are used
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.02 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.02 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [5 5]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
par.ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'median'; 'ssim-ml';'median';'entropy';'none';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
% Correlation of projections and flat fields. Essential for DCM data. Typically improves reconstruction quality of DMM data, too.
% Available methods ('ssim-ml'/'entropy' usually work best):
% 'none' : no correlation, for DPC
% 'mean'/'median': mean/median flat
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
image_correlation.num_flats = 11; % integer. number of best maching flat fields used for correction
image_correlation.area_width = [1 100];% 2-vector. correlation area: index vector or relative/absolute position of [first pix, last pix], negative indexing is supported
image_correlation.area_height = [0.25 0.75]; % 2/vector. correlation area [bottom top]: index vector or relative/absolute position of [first pix, last pix]
image_correlation.filter = 1; % bool, filter ROI before correlation
image_correlation.filter_type = 'median'; % string. correlation ROI filter type, currently only 'median' is implemnted
image_correlation.filter_parameter = {[5 5], 'symmetric'}; % cell. filter paramaters to be parsed with {:}
% 'median' : using medfilt2, parameters: {[M N]-neighboorhood, 'padding'}
% 'wiener' : using wiender2, parameters: {[M N]-neighboorhood}
ring_filter.apply = 1; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft_dec_levels = 1:6; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft_wname = 'db7';'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft_sigma = 3; % integer scalar. suppression factor for 'wavelet-fft'
ring_filter.jm_median_width = 11; % integer scalar or vector. median averaging filter to be applied to angular averaged sinogram, multiple widths are applied consecutively, eg [3 11 21 31 39];
par.strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
par.norm_sino = 0; % not recommended, can introduce severe artifacts, but sometimes improves quality
% Workaround correction for image distortions using a quadratic dilation/compression of the projections/sinogram
% Preferably, projection cropping of laterally shifted projection 'crop_proj' is not used
par.distortion_correction_distance = 0; % scalar, in binned pixel, distance between two regions in the tomogram that can be properly reconstructed using different rotation axis offsets, if 0: no correction done
par.distortion_correction_outer_offset = 0; % scalar, in pixel, rotation axis offset for the outer region. the offset for the inner region is used for reconstruction
par.distortion_correction_exponent = 2; % scalar,  exponent of interpolation function: xq = x - 2 * offset_diff * (x / dist_offset).^exponent;
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % bool. See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % bool. before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % bool. Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie'; % string. Available methods: 'qp' 'ctf' 'tie' 'qp2' 'qpcut', 'tieNLO_Schwinger'
% Interactive phase retrieval not supported for method 'tieNLO_Schwinger'
phase_retrieval.reg_par = 1.0; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
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
tomo.reco_mode = '3D';'slice'; % string. slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = [];%[-1 1 -1 1 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax] for excentric rot axis pos or extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
% Orientation using Matlab's matrix notation: relative coordinates in a horizontal reconstruction plane using imagej (-0.5,-0.5) = top left
% pixel, (0.5,0.5) = bottom right, (0,0) = center
% [left, right, top, bottom, bottom slice, top slice]
tomo.vol_shape = []; %[1 1 1] integer vector. shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_full_range = [];% 2 * pi ;[]; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle_offset = 0; % global rotation of reconstructed volume
tomo.interpolate_missing_angles = 0; % limited or missing angle tomography
tomo.rot_axis_offset = [] / 1 * par.raw_bin; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.vert_shift = []; % vertical shift for spiral/helical CT
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
tomo.rot_axis_search_auto = 0; % find extrema of metric within search range
tomo.rot_axis_search_range = []; % search reach for automatic determination of the rotation axis offset, overwrite interactive result if not empty
tomo.rot_axis_search_metric = 'iso-grad'; % string: 'neg','entropy','iso-grad','laplacian','entropy-ML','abs'. Metric to find rotation axis offset
tomo.rot_axis_search_extrema = 'max'; % string: 'min'/'max'. chose min or maximum position
tomo.rot_axis_search_fit = 1; % bool: fit calculated metrics and find extrema, otherwise use extrema from search range
tomo.rot_axis_offset_metric_roi = []; % 4-vector: [. ROI for metric calculation. roi = [y0, x0, y1-y0, x1-x0]. (x,y)=(0,0)=upper left
tomo.rot_axis_search_slice = []; % scalar: slice used to find rot axis. if empty: uses slice from interactive mode, if that is empty uses central slice.
tomo.rot_axis_search_range_from_interactive = 0; % boolean: use search range from interactive mode
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % string. absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID, mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = ''; % parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder_flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder_phase_map = ''; % subfolder in 'phase_map'
write.subfolder_sino = ''; % subfolder in 'sino'
write.subfolder_reco = ''; % subfolder in 'reco'
write.flatcor = 1; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.float_adapthisteq = 0; % save float with adaptive histogram equalization filter. only for 3D recos currently
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression_method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression_method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 0; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression_method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression_method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression_method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compress dynamic range into [0, 1]
write.compression_parameter = [0.02 0.02]; % compression-method specific parameters
% methods for the compression of the dynamic range:
% 'outlier' : [LOW, HIGH] = write.compression_parameter, eg. [0.01 0.03], outlier given in percent, if scalar LOW = HIGH.
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression_parameter, eg. [-0.01 1]
% 'std' : NUM = write.compression_parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression_parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__
write.outputformat = 'tif';'hdf_volume'; % string. Not yet implemented for all reco modes
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 1; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 0; % find laminography tilt instead camera tilt
interactive_mode.angles = 0; % reconstruct slices with different scalings of angles
interactive_mode.angle_scaling_default_search_range = []; % if empty: use a variaton of -/+5 * (angle increment / maximum angle)
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 1; % Interactive retrieval to determine regularization parameter
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
interactive_mode.show_stack_imagej = 1; % use imagej instead of MATLAB to scroll through images during interactive mode
interactive_mode.show_stack_imagej_use_virtual = 1; % use virtual stack for faster loading, but slower scrolling
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.astra_link_data = 1; % boolean: ASTRA data objects become references to Matlab arrays. Reduces memory issues.
par.gpu_index = []; % integer vector: indices of GPU devices to use, Matlab notation: index starts from 1. default: [], uses all
par.use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.use_gpu_in_parfor = 1; % boolean
pixel_filter_sino.use_gpu = par.use_gpu_in_parfor;
par.poolsize = 0.5; % scalar: number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.5; % scalar: elative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
phase_retrieval.use_parpool = 1; % bool. Disable parpool when out-of-memory error occurs during phase retrieval.
par.window_state = 'minimized';'normal';'maximized';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proc_path = '/asap3/petra3/gpfs/p07/2023/data/11017598/processed/';
raw_path =  '/asap3/petra3/gpfs/p07/2023/data/11017598/raw/';

par.scan_path = [proc_path 'tig002_81217_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
ADD

par.scan_path = [proc_path 'tig003_84785_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
par.read_sino = 1; 
par.read_sino_folder = 'trans02_180';
par.read_sino_trafo = @(x) (x);
par.sino_roi = []; 
par.filter_sino = 1;
pixel_filter_sino.threshold_hot = 0;0.0005;
pixel_filter_sino.threshold_dark = 0;0.00005;
pixel_filter_sino.medfilt_neighboorhood = [3 3];
pixel_filter_sino.filter_dead_pixel = 1;
pixel_filter_sino.filter_Inf = 1;
pixel_filter_sino.filter_NaN = 1;
pixel_filter_sino.verbose = 0;
pixel_filter_sino.use_gpu = par.use_gpu_in_parfor;
ring_filter.apply = 1;
tomo.rot_axis_offset = 0;
interactive_mode.rot_axis_pos = 0;
phase_retrieval.apply = 1;
phase_retrieval.reg_par = 1.0;
interactive_mode.phase_retrieval = 0;
ADD

par.scan_path = [proc_path 'tig004_81166_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig005_84243_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig006_81217_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig008_81162_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig009_84244_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig010_80513_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig011_80510_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig012_84786_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig013_81954_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig014_84783_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig015_80610_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig016_84787_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig017_85680_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

par.scan_path = [proc_path 'tig018_81951_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

% Binning 1
par.scan_path = [proc_path 'tig004_81166_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
par.raw_bin = 1;
par.read_sino_folder = 'trans01_180';
%fprintf('\n%s: %u',par.nexus_path, exist(par.nexus_path,'dir'))
phase_retrieval.apply = 0;
ADD
phase_retrieval.apply = 1;
ADD

% Different phase reg par
par.scan_path = [proc_path 'tig004_81166_mousebrain_202309_DESY'];
par.nexus_path = [regexprep(par.scan_path,'processed','raw') '_height_a'];
par.raw_bin = 2;
par.read_sino_folder = 'trans02_180';
phase_retrieval.apply = 1;
phase_retrieval.reg_par = 0.0;ADD
phase_retrieval.reg_par = 0.2;ADD
phase_retrieval.reg_par = 0.4;ADD
phase_retrieval.reg_par = 0.6;ADD
phase_retrieval.reg_par = 0.8;ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
