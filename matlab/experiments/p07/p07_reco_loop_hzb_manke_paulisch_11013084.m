function p07_reco_loop_hzb_manke_paulisch_11013084( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
% Created on 24-May-2022 by moosmanj

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


pp_parameter_switch % DO NOT DELETE THIS LINE

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%
par.scan_path = pwd; % string/pwd. pwd: change to directory of the scan to be reconstructed, string: absolute scan path, last_folder_modified('folder')
par.ref_path = {}; % cell of strings. Additonal data sets to be included for the correlation of projections and reference images
par.read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
par.read_flatcor_path = ''; % absolute path containing flat-field corrected projections
par.read_flatcor_trafo = @(im) im; %fliplr( im ); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
par.read_sino_folder = ''; % subfolder to scan path
par.read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.energy = []; % eV! if empty: read from log file (log file values can be ambiguous or even missing sometimes)
par.sample_detector_distance = []; % in m. if empty: read from log file
par.eff_pixel_size = []; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
par.pixel_scaling = []; % to account for mismatch of eff_pixel_size with, ONLY APPLIED BEFORE TOMOGRAPHIC RECONSTRUCTION, HAS TO BE CHANGED!
par.read_image_log = 0; % bool, default: 0. Read metadata from image log instead hdf5, if image log exists
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.raw_bin = 3; % projection binning factor: integer
par.raw_roi = []; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
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
par.ref_range = [];%1:100;[]; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
par.wo_crop = 0; % Do not crop images to account for random lateral shift
par.virt_s_pos = 0; % Correct sample position in reconsructed volume if virtual sample position motors are used
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.02 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.02 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [5 5]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
par.ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'median';'entropy';'none';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
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
image_correlation.num_flats = 11; % number of best maching flat fields used for correction
image_correlation.area_width = [1 100];%[-100 1];% correlation area: index vector or relative/absolute position of [first pix, last pix], negative indexing is supported
image_correlation.area_height = [0.6 1.0]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
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
phase_retrieval.reg_par = 2; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
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
tomo.vol_size = [];[-1 1 -1 1 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax]for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
% [ slice left,  slice right,  slice top,  slice bottom, bottom slice, top slice], 2D Matlab notation: (0,0)=topleft:
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_full_range = [];% 2 * pi ;[]; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle_offset = 0; % global rotation of reconstructed volume
tomo.interpolate_missing_angles = 0; % limited or missing angle tomography
tomo.rot_axis_offset = [] / 2 * par.raw_bin; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR posITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.vert_shift = [];
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
tomo.rot_axis_search_metric = 'neg'; % string: 'neg','entropy','iso-grad','laplacian','entropy-ML','abs'. Metric to find rotation axis offset
tomo.rot_axis_search_extrema = 'max'; % string: 'min'/'max'. chose min or maximum position
tomo.rot_axis_search_fit = 1; % bool: fit calculated metrics and find extrema, otherwise use extrema from search range
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
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
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression_method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression_method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 1; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression_method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression_method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression_method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression_parameter = [0.02 0.02]; % compression-method specific parameters
% dynamic range is compressed s.t. new dynamic range assumes
% 'outlier' : [LOW, HIGH] = write.compression_parameter, eg. [0.01 0.03], outlier given in percent, if scalear LOW = HIGH.
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression_parameter, eg. [-0.01 1]
% 'std' : NUM = write.compression_parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression_parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 1; % show images and plots during reconstruction
par.skip_gpu_info = 1;
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
par.poolsize = 0.8; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.8; % Relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
par.gpu_index = tomo.astra_gpu_index;
par.use_gpu_in_parfor = 1;
phase_retrieval.use_parpool = 1; % bool. Disable parpool when out-of-memory error occurs during phase retrieval.
par.window_state = 'minimized';'normal';'maximized';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p07/2022/data/11013084/raw/';


par.raw_bin = 3;
par.raw_roi = -2;
ring_filter.apply = 0;
phase_retrieval.apply = 0;
tomo.rot_axis_search_auto = 1;
tomo.rot_axis_search_range = -12:4;
tomo.rot_axis_search_metric = 'entropy';
tomo.rot_axis_search_extrema = 'min';
tomo.rot_axis_search_fit = 1;
tomo.slice = 0.4;
write.to_scratch = 0;
interactive_mode.rot_axis_pos = 0;
interactive_mode.slice_number = 0.4;
% END OF QUICK SWITCH TO ALTERNATIVE SET OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

par.scan_path = [raw_path 'for001_ag97_trocken']; ADD

par.scan_path = [raw_path 'for002_ag97_127_4a']; ADD

par.scan_path = [raw_path 'for003_ag97_127_4a_step0']; ADD

par.scan_path = [raw_path 'for004_ag97_135_1a_step0']; ADD

par.scan_path = [raw_path 'for005_ag97_135_1a_step1']; ADD

par.scan_path = [raw_path 'for006_ag97_135_1b_step0']; ADD

par.scan_path = [raw_path 'for007_ag97_135_1c_step0']; ADD

par.scan_path = [raw_path 'for008_ag97_135_1c_step1']; ADD

par.scan_path = [raw_path 'for009_ag97_135_1c_step2']; ADD

par.scan_path = [raw_path 'for010_ag97_135_1d_step0']; ADD

par.scan_path = [raw_path 'for011_ag97_135_1d_step1']; ADD

par.scan_path = [raw_path 'for012_ag97_135_1d_step2']; ADD

par.scan_path = [raw_path 'for013_ag97_135_1d_step3']; ADD

par.scan_path = [raw_path 'for014_ag97_135_1d_step4']; ADD

par.scan_path = [raw_path 'for015_ag97_135_1e_step0']; ADD

par.scan_path = [raw_path 'for016_ag97_135_1e_step1']; ADD

par.scan_path = [raw_path 'for017_ag97_135_1f_step0']; ADD

par.scan_path = [raw_path 'for018_ag97_135_1g_step0']; ADD

par.scan_path = [raw_path 'for019_ag97_135_1g1_step0']; ADD

par.scan_path = [raw_path 'for020_ag97_135_1h_step0']; ADD

par.scan_path = [raw_path 'for021_ag97_135_1h_step1']; ADD

par.scan_path = [raw_path 'for022_ag97_135_1h_step2']; ADD

par.scan_path = [raw_path 'for023_ag97_135_1h_step3']; ADD

par.scan_path = [raw_path 'for024_ag97_135_1i_step0']; ADD

par.scan_path = [raw_path 'for025_ag97_135_1j_step0']; ADD

par.scan_path = [raw_path 'for026_ag97_135_1j_step1']; ADD

par.scan_path = [raw_path 'for027_ag97_135_1j_step2']; ADD

par.scan_path = [raw_path 'for028_ag97_135_1k_step0']; ADD

par.scan_path = [raw_path 'for029_ag97_135_1k_step1']; ADD

par.scan_path = [raw_path 'for030_ag97_135_1k_step2']; ADD

par.scan_path = [raw_path 'for031_ag97_135_1k_step3']; ADD

par.scan_path = [raw_path 'for032_ag97_135_1l_step0']; ADD

par.scan_path = [raw_path 'for033_ag97_135_1l_step1']; ADD

par.scan_path = [raw_path 'for034_ag97_135_1l_step2']; ADD

par.scan_path = [raw_path 'for035_ag97_135_1l_step3']; ADD

par.scan_path = [raw_path 'for036_ag97_135_1l_step4']; ADD

par.scan_path = [raw_path 'for037_ag97_135_1l_step5']; ADD

par.scan_path = [raw_path 'for038_ag97_135_1l_step6']; ADD

par.scan_path = [raw_path 'for039_ag97_135_1l_step7']; ADD

par.scan_path = [raw_path 'for040_ag97_135_1l_step8']; ADD

par.scan_path = [raw_path 'for041_ag97_135_1l_step9']; ADD

par.scan_path = [raw_path 'for042_ag97_135_1m_step0']; ADD

par.scan_path = [raw_path 'for043_ag97_135_1m_step01']; ADD

par.scan_path = [raw_path 'for044_ag97_135_1m_step01']; ADD

par.scan_path = [raw_path 'for045_ag97_135_1m_step02']; ADD

par.scan_path = [raw_path 'for046_ag97_135_1m_step03']; ADD

par.scan_path = [raw_path 'for047_ag97_135_1m_step04']; ADD

par.scan_path = [raw_path 'for048_ag97_135_1m_step05']; ADD

par.scan_path = [raw_path 'for049_ag97_135_1m_step06']; ADD

par.scan_path = [raw_path 'for050_ag97_135_1m_step07']; ADD

par.scan_path = [raw_path 'for051_ag97_135_1m_step08']; ADD

par.scan_path = [raw_path 'for052_ag97_135_1m_step09']; ADD

par.scan_path = [raw_path 'for053_ag97_135_1m_step10']; ADD

par.scan_path = [raw_path 'for054_ag97_135_1n_step0']; ADD

par.scan_path = [raw_path 'for055_ag97_135_1n_step01']; ADD

par.scan_path = [raw_path 'for056_ag97_135_1n_step02']; ADD

par.scan_path = [raw_path 'for057_ag97_135_1n_step03']; ADD

par.scan_path = [raw_path 'for058_ag97_135_1n_step04']; ADD

par.scan_path = [raw_path 'for059_ag97_135_1n_step05']; ADD

par.scan_path = [raw_path 'for060_ag92_141_2a_step0']; ADD

par.scan_path = [raw_path 'for061_ag92_141_2a_step1']; ADD

par.scan_path = [raw_path 'for062_ag92_141_2a_step2']; ADD

par.scan_path = [raw_path 'for063_ag92_141_2a_step3']; ADD

par.scan_path = [raw_path 'for064_ag92_141_2a_step4']; ADD

par.scan_path = [raw_path 'for065_ag92_141_2a_step5']; ADD

par.scan_path = [raw_path 'for066_ag92_141_2a_step6']; ADD

par.scan_path = [raw_path 'for067_ag92_141_2a_step7']; ADD

par.scan_path = [raw_path 'for068_ag92_141_2a_step8']; ADD

par.scan_path = [raw_path 'for069_ag92_141_2a_step9']; ADD

par.scan_path = [raw_path 'for070_ag92_141_2b_step0']; ADD

par.scan_path = [raw_path 'for071_ag97_141_2b_step01']; ADD

par.scan_path = [raw_path 'for072_ag97_141_2b_step02']; ADD

par.scan_path = [raw_path 'for073_ag97_141_2b_step03']; ADD

par.scan_path = [raw_path 'for074_ag97_141_2b_step04']; ADD

par.scan_path = [raw_path 'for075_ag97_141_2b_step05']; ADD

par.scan_path = [raw_path 'for076_ag97_135_1o_step0']; ADD

par.scan_path = [raw_path 'for077_ag97_135_1o_step1']; ADD

par.scan_path = [raw_path 'for078_ag97_135_1o_step2']; ADD

par.scan_path = [raw_path 'for079_ag97_135_1o_step3']; ADD

par.scan_path = [raw_path 'for080_ag97_135_1o_step4']; ADD

par.scan_path = [raw_path 'for081_ag97_135_1o_step5']; ADD

par.scan_path = [raw_path 'for082_ag97_135_1o_step6']; ADD

par.scan_path = [raw_path 'for083_ag97_135_1o_step7']; ADD

par.scan_path = [raw_path 'for084_ag97_135_1o_step8']; ADD

par.scan_path = [raw_path 'for085_ag97_135_1o_step9']; ADD

par.scan_path = [raw_path 'for086_ag97_135_1o_step10']; ADD

par.scan_path = [raw_path 'for087_ag97_135_1o_step11']; ADD

par.scan_path = [raw_path 'for088_ag92_141_2c_step0']; ADD

par.scan_path = [raw_path 'for089_ag92_141_2c_step1']; ADD

par.scan_path = [raw_path 'for090_ag92_141_2c_step2']; ADD

par.scan_path = [raw_path 'for091_ag92_141_2c_step3']; ADD

par.scan_path = [raw_path 'for092_ag92_141_2c_step4']; ADD

par.scan_path = [raw_path 'for093_ag92_141_2c_step5']; ADD

par.scan_path = [raw_path 'for094_ag92_141_2c_step6']; ADD

par.scan_path = [raw_path 'for095_ag92_141_2d_step0']; ADD

par.scan_path = [raw_path 'for096_ag92_141_2d_step01']; ADD

par.scan_path = [raw_path 'for097_ag92_141_2d_step02']; ADD

par.scan_path = [raw_path 'for098_ag92_141_2d_step03']; ADD

par.scan_path = [raw_path 'for099_ag92_141_2d_step04']; ADD

par.scan_path = [raw_path 'for100_ag92_141_2d_step05']; ADD

par.scan_path = [raw_path 'for101_ag92_141_2e_step0']; ADD

par.scan_path = [raw_path 'for102_ag92_141_2e_step1']; ADD

par.scan_path = [raw_path 'for103_ag92_141_2e_step2']; ADD

par.scan_path = [raw_path 'for104_ag92_141_2e_step3']; ADD

par.scan_path = [raw_path 'for105_ag92_141_2e_step4']; ADD

par.scan_path = [raw_path 'for106_ag92_141_2e_step5']; ADD

par.scan_path = [raw_path 'for107_ag92_141_2e_step6']; ADD

par.scan_path = [raw_path 'for108_ag92_141_2f_step0']; ADD

par.scan_path = [raw_path 'for109_ag92_141_2f_step01']; ADD

par.scan_path = [raw_path 'for110_ag92_141_2f_step02']; ADD

par.scan_path = [raw_path 'for111_ag92_141_2f_step03']; ADD

par.scan_path = [raw_path 'for112_ag92_141_2f_step04']; ADD

par.scan_path = [raw_path 'for113_ag92_141_2f_step05']; ADD

par.scan_path = [raw_path 'for113a_test für reco']; ADD

par.scan_path = [raw_path 'for113a_test für reco (multy_2)']; ADD

par.scan_path = [raw_path 'for114_ag95_138_3a_step0']; ADD

par.scan_path = [raw_path 'for115_ag95_138_3a_step1']; ADD

par.scan_path = [raw_path 'for116_ag95_138_3a_step2']; ADD

par.scan_path = [raw_path 'for117_ag95_138_3a_step2']; ADD

par.scan_path = [raw_path 'for118_ag95_138_3a_step2']; ADD

par.scan_path = [raw_path 'for119_ag95_138_3a_step5']; ADD

par.scan_path = [raw_path 'for120_ag95_138_3a_step6']; ADD

par.scan_path = [raw_path 'for121_ag95_138_3a_step7']; ADD

par.scan_path = [raw_path 'for122_ag95_138_3a_step8']; ADD

par.scan_path = [raw_path 'for123_ag95_138_3a_step9']; ADD

par.scan_path = [raw_path 'for124_ag95_138_3a_step10']; ADD

par.scan_path = [raw_path 'for125_ag95_138_3a_step11']; ADD

par.scan_path = [raw_path 'for126_ag95_138_3b_step0']; ADD

par.scan_path = [raw_path 'for127_ag95_138_3b_step01']; ADD

par.scan_path = [raw_path 'for128_ag95_138_3b_step02']; ADD

par.scan_path = [raw_path 'for129_ag95_138_3b_step03']; ADD

par.scan_path = [raw_path 'for130_ag95_138_3b_step04']; ADD

par.scan_path = [raw_path 'for131_ag95_138_3b_step05']; ADD

par.scan_path = [raw_path 'for132_ag95_138_3c_step0']; ADD

par.scan_path = [raw_path 'for133_ag95_138_3c_step1']; ADD

par.scan_path = [raw_path 'for134_ag95_138_3c_step2']; ADD

par.scan_path = [raw_path 'for135_ag95_138_3c_step3']; ADD

par.scan_path = [raw_path 'for136_ag95_138_3c_step4']; ADD

par.scan_path = [raw_path 'for137_ag95_138_3c_step5']; ADD

par.scan_path = [raw_path 'for138_ag95_138_3c_step6']; ADD

par.scan_path = [raw_path 'for139_ag95_138_3c_step7']; ADD

par.scan_path = [raw_path 'for140_ag95_138_3c_step8']; ADD

par.scan_path = [raw_path 'for141_ag95_138_3c_step9']; ADD

par.scan_path = [raw_path 'for142_ag95_138_3d_step0']; ADD

par.scan_path = [raw_path 'for143_ag95_138_3d_step01']; ADD

par.scan_path = [raw_path 'for144_ag95_138_3d_step02']; ADD

par.scan_path = [raw_path 'for145_ag95_138_3d_step03']; ADD

par.scan_path = [raw_path 'for146_ag95_138_3d_step04']; ADD

par.scan_path = [raw_path 'for147_ag95_138_3d_step05']; ADD

par.scan_path = [raw_path 'for148_ag95_138_3e_step0']; ADD

par.scan_path = [raw_path 'for149_ag95_138_3e_step1']; ADD

par.scan_path = [raw_path 'for150_ag95_138_3e_step2']; ADD

par.scan_path = [raw_path 'for151_ag95_138_3e_step3']; ADD

par.scan_path = [raw_path 'for152_ag95_138_3e_step4']; ADD

par.scan_path = [raw_path 'for153_ag95_138_3e_step5']; ADD

par.scan_path = [raw_path 'for154_ag95_138_3e_step6']; ADD

par.scan_path = [raw_path 'for155_ag95_138_3f_step0']; ADD

par.scan_path = [raw_path 'for156_ag95_138_3f_step01']; ADD

par.scan_path = [raw_path 'for157_ag95_138_3f_step02']; ADD

par.scan_path = [raw_path 'for158_ag95_138_3f_step03']; ADD

par.scan_path = [raw_path 'for159_ag95_138_3f_step04']; ADD

par.scan_path = [raw_path 'for160_ag95_138_3f_step05']; ADD

par.scan_path = [raw_path 'for161_ag97_135_1p_step0']; ADD

par.scan_path = [raw_path 'for162_ag97_135_1p_step01']; ADD

par.scan_path = [raw_path 'for163_ag97_135_1p_step02']; ADD

par.scan_path = [raw_path 'for164_ag97_135_1p_step03']; ADD

par.scan_path = [raw_path 'for165_ag97_135_1p_step04']; ADD

par.scan_path = [raw_path 'for166_ag97_135_1p_step05']; ADD

par.scan_path = [raw_path 'for167_ag97_135_1q_step0']; ADD

par.scan_path = [raw_path 'for168_ag97_135_1q_step1']; ADD

par.scan_path = [raw_path 'for169_ag97_135_1q_step2']; ADD

par.scan_path = [raw_path 'for170_ag97_135_1q_step3']; ADD

par.scan_path = [raw_path 'for171_ag97_135_1q_step4']; ADD

par.scan_path = [raw_path 'for172_ag97_135_1q_step5']; ADD

par.scan_path = [raw_path 'for173_ag97_135_1q_step6']; ADD

par.scan_path = [raw_path 'for174_ag97_135_1r_step0']; ADD

par.scan_path = [raw_path 'for175_ag97_135_1r_step1']; ADD

par.scan_path = [raw_path 'for176_ag97_135_1r_step2']; ADD

par.scan_path = [raw_path 'for177_ag97_135_1r_step3']; ADD

par.scan_path = [raw_path 'for178_ag97_135_1r_step4']; ADD

par.scan_path = [raw_path 'for179_ag97_135_1r_step5']; ADD

par.scan_path = [raw_path 'for180_ag97_135_1r_step6']; ADD

par.scan_path = [raw_path 'for181_ag97ss_140_9a_step0']; ADD

par.scan_path = [raw_path 'for182_ag97ss_140_9a_step1']; ADD

par.scan_path = [raw_path 'for183_ag97ss_140_9a_step2']; ADD

par.scan_path = [raw_path 'for184_ag97ss_140_9a_step3']; ADD

par.scan_path = [raw_path 'for185_ag97ss_140_9a_step4']; ADD

par.scan_path = [raw_path 'for186_ag97ss_140_9a_step5']; ADD

par.scan_path = [raw_path 'for187_ag97ss_140_9a_step6']; ADD

par.scan_path = [raw_path 'for188_ag97_167_1a_co2rr_step0']; ADD

par.scan_path = [raw_path 'for189_ag97_167_1a_co2rr_step1']; ADD

par.scan_path = [raw_path 'for190_ag97_167_1a_co2rr_step2']; ADD

par.scan_path = [raw_path 'for191_ag97_167_1a_co2rr_step3']; ADD

par.scan_path = [raw_path 'for192_ag97_167_1a_co2rr_step4']; ADD

par.scan_path = [raw_path 'for193_ag97_167_1a_co2rr_step5']; ADD

par.scan_path = [raw_path 'for194_ag97_167_1a_co2rr_step6']; ADD

par.scan_path = [raw_path 'for195_ag97_167_1b_co2rr_step0']; ADD

par.scan_path = [raw_path 'for196_ag97_167_1b_step01']; ADD

par.scan_path = [raw_path 'for197_ag97_167_1b_step02']; ADD

par.scan_path = [raw_path 'for198_ag97_167_1b_step03']; ADD

par.scan_path = [raw_path 'for199_ag97_167_1b_step04']; ADD

par.scan_path = [raw_path 'for200_ag97_167_1b_step05']; ADD

par.scan_path = [raw_path 'for201_ag97_167_1c_co2rr_step0']; ADD

par.scan_path = [raw_path 'for202_ag97_167_1c_co2rr_step1']; ADD

par.scan_path = [raw_path 'for203_ag97_167_1c_co2rr_step2']; ADD

par.scan_path = [raw_path 'for204_ag97_167_1c_co2rr_step3']; ADD

par.scan_path = [raw_path 'for205_ag97_167_1c_co2rr_step4']; ADD

par.scan_path = [raw_path 'for206_ag97_167_1c_co2rr_step5']; ADD

par.scan_path = [raw_path 'for207_ag97_167_1c_co2rr_step6']; ADD

par.scan_path = [raw_path 'for208_ag97_167_1d_co2rr_step0']; ADD

par.scan_path = [raw_path 'for209_ag97_167_1d_step01']; ADD

par.scan_path = [raw_path 'for210_ag97_167_1d_step02']; ADD

par.scan_path = [raw_path 'for211_ag97_167_1d_step03']; ADD

par.scan_path = [raw_path 'for212_ag97_167_1d_step04']; ADD

par.scan_path = [raw_path 'for213_ag97_167_1d_step05']; ADD

par.scan_path = [raw_path 'for214_ag97_167_1e_co2rr_step0']; ADD

par.scan_path = [raw_path 'for215_ag97_167_1e_co2rr_step1']; ADD

par.scan_path = [raw_path 'for216_ag97_167_1e_co2rr_step2']; ADD

par.scan_path = [raw_path 'for217_ag97_167_1e_co2rr_step3']; ADD

par.scan_path = [raw_path 'for218_ag97_167_1e_co2rr_step4']; ADD

par.scan_path = [raw_path 'for219_ag97_167_1e_co2rr_step5']; ADD

par.scan_path = [raw_path 'for220_ag97_167_1e_co2rr_step6']; ADD

par.scan_path = [raw_path 'for221_ag97_167_1f_co2rr_step0']; ADD

par.scan_path = [raw_path 'for222_ag97_167_1f_step01']; ADD

par.scan_path = [raw_path 'for223_ag97_167_1f_step02']; ADD

par.scan_path = [raw_path 'for224_ag97_167_1f_step03']; ADD

par.scan_path = [raw_path 'for225_ag97_167_1f_step04']; ADD

par.scan_path = [raw_path 'for226_ag97_167_1f_step05']; ADD

par.scan_path = [raw_path 'for227_ag92_141_2g_co2rr_step0']; ADD

par.scan_path = [raw_path 'for228_ag92_141_2g_co2rr_step01']; ADD

par.scan_path = [raw_path 'for229_ag92_141_2g_co2rr_step02']; ADD

par.scan_path = [raw_path 'for230_ag92_141_2g_co2rr_step03']; ADD

par.scan_path = [raw_path 'for231_ag92_141_2g_co2rr_step04']; ADD

par.scan_path = [raw_path 'for232_ag92_141_2g_co2rr_step05']; ADD

par.scan_path = [raw_path 'for233_ag92_141_2h_co2rr_step0']; ADD

par.scan_path = [raw_path 'for234_ag92_141_2h_co2rr_step1']; ADD

par.scan_path = [raw_path 'for235_ag92_141_2h_co2rr_step2']; ADD

par.scan_path = [raw_path 'for236_ag92_141_2h_co2rr_step3']; ADD

par.scan_path = [raw_path 'for237_ag92_141_2h_co2rr_step4']; ADD

par.scan_path = [raw_path 'for238_ag92_141_2h_co2rr_step5']; ADD

par.scan_path = [raw_path 'for239_ag92_141_2h_co2rr_step6']; ADD

par.scan_path = [raw_path 'for240_ag92_141_2i_co2rr_step0']; ADD

par.scan_path = [raw_path 'for241_ag92_141_2i_co2rr_step1']; ADD

par.scan_path = [raw_path 'for242_ag92_141_2i_co2rr_step2']; ADD

par.scan_path = [raw_path 'for243_ag92_141_2i_co2rr_step']; ADD

par.scan_path = [raw_path 'for244_ag92_141_2i_co2rr_step4']; ADD

par.scan_path = [raw_path 'for245_ag92_141_2i_co2rr_step5']; ADD

par.scan_path = [raw_path 'for246_ag92_141_2i_co2rr_step6']; ADD

par.scan_path = [raw_path 'for247_ag92_141_2j_co2rr_step0']; ADD

par.scan_path = [raw_path 'for248_ag92_141_2j_co2rr_step01']; ADD

par.scan_path = [raw_path 'for249_ag92_141_2j_co2rr_step02']; ADD

par.scan_path = [raw_path 'for250_ag92_141_2j_co2rr_step03']; ADD

par.scan_path = [raw_path 'for251_ag92_141_2j_co2rr_step04']; ADD

par.scan_path = [raw_path 'for252_ag92_141_2j_co2rr_step05']; ADD

par.scan_path = [raw_path 'for253_ag95_138_3g_co2rr_step0']; ADD

par.scan_path = [raw_path 'for254_ag95_138_3g_co2rr_step1']; ADD

par.scan_path = [raw_path 'for255_ag95_138_3g_co2rr_step2']; ADD

par.scan_path = [raw_path 'for256_ag95_138_3g_co2rr_step3']; ADD

par.scan_path = [raw_path 'for257_ag95_138_3g_co2rr_step4']; ADD

par.scan_path = [raw_path 'for258_ag95_138_3g_co2rr_step5']; ADD

par.scan_path = [raw_path 'for259_ag95_138_3g_co2rr_step6']; ADD

par.scan_path = [raw_path 'for260_ag95_138_3h_co2rr_step0']; ADD

par.scan_path = [raw_path 'for261_ag95_138_3h_co2rr_step01']; ADD

par.scan_path = [raw_path 'for262_ag95_138_3h_co2rr_step02']; ADD

par.scan_path = [raw_path 'for263_ag95_138_3h_co2rr_step03']; ADD

par.scan_path = [raw_path 'for264_ag95_138_3h_co2rr_step04']; ADD

par.scan_path = [raw_path 'for265_ag95_138_3h_co2rr_step05']; ADD

par.scan_path = [raw_path 'for266_ag95_138_3i_co2rr_step0']; ADD

par.scan_path = [raw_path 'for267_ag95_138_3i_co2rr_step01']; ADD

par.scan_path = [raw_path 'for268_ag95_138_3i_co2rr_step02']; ADD

par.scan_path = [raw_path 'for269_ag95_138_3i_co2rr_step03']; ADD

par.scan_path = [raw_path 'for270_ag95_138_3i_co2rr_step04']; ADD

par.scan_path = [raw_path 'for271_ag95_138_3i_co2rr_step05']; ADD

par.scan_path = [raw_path 'for272_ag95_138_3j_co2rr_step0']; ADD

par.scan_path = [raw_path 'for273_ag95_138_3j_co2rr_step1']; ADD

par.scan_path = [raw_path 'for274_ag95_138_3j_co2rr_step2']; ADD

par.scan_path = [raw_path 'for275_ag95_138_3j_co2rr_step3']; ADD

par.scan_path = [raw_path 'for276_ag95_138_3j_co2rr_step4']; ADD

par.scan_path = [raw_path 'for277_ag95_138_3j_co2rr_step5']; ADD

par.scan_path = [raw_path 'for278_ag95_138_3j_co2rr_step6']; ADD

par.scan_path = [raw_path 'for279_ag92_141_2k_co2rr_step0']; ADD

par.scan_path = [raw_path 'for280_ag92_141_2k_co2rr_step01']; ADD

par.scan_path = [raw_path 'for281_ag92_141_2k_co2rr_step02']; ADD

par.scan_path = [raw_path 'for282_ag92_141_2k_co2rr_step03']; ADD

par.scan_path = [raw_path 'for283_ag92_141_2k_co2rr_step04']; ADD

par.scan_path = [raw_path 'for284_ag92_141_2k_co2rr_step05']; ADD

par.scan_path = [raw_path 'for285_kuevette_1p5M_CsHCO3']; ADD

par.scan_path = [raw_path 'for286_kuevette_1p5M_CsHCO3']; ADD

par.scan_path = [raw_path 'for287_kuevette_50v50_50wtproz_CsOH_30wtproz_NaOH']; ADD

par.scan_path = [raw_path 'for288_kuevette_50v50_50wtproz_CsOH_30wtproz_NaOH']; ADD

par.scan_path = [raw_path 'for289_cov_a_step0']; ADD

par.scan_path = [raw_path 'for290_cov_a_step1']; ADD

par.scan_path = [raw_path 'for291_cov_a_step2']; ADD

par.scan_path = [raw_path 'for292_cov_a_step3']; ADD

par.scan_path = [raw_path 'for293_cov_a_step4']; ADD

par.scan_path = [raw_path 'for294_cov_a_step5']; ADD

par.scan_path = [raw_path 'for295_cov_a_step6']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
