function p05_reco_loop_synchroload2017oct_11003440_001( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
% PARAMETER section with 'SET_DEFAULT', then the first time 'ADD' is called
% defines default parameters.
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
% Created on 19-Mar-2018 by moosmanj

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

read_flatcor = 0; % read flatfield-corrected images from disc, skips preprocessing
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = [];
raw_bin = 2; 
bin_before_filtering = 0;
im_trafo = ''; 
excentric_rot_axis = 0;
crop_at_rot_axis = 0;
stitch_projections = 0;
stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
proj_range = 1; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = 1; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixe
ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';
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
strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
decimal_round_precision = 2; % precision when rounding pixel shifts
% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0;
phase_retrieval.apply_before = 0;
phase_retrieval.post_binning_factor = 1;
phase_retrieval.method = 'tie';'qpcut';
phase_retrieval.reg_par = 1.5;
phase_retrieval.bin_filt = 0.15;
phase_retrieval.cutoff_frequ = 2 * pi;
phase_retrieval.padding = 1;
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, intendet to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D'; 'slice'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = []; %[-0.5 0.5 -0.5 0.5 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle.full_range = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle.offset = pi; % global rotation of reconstructed volume
tomo.rot_axis.offset = 0;
tomo.rot_axis.offset_shift_range = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.rot_axis.position = [];
tomo.rot_axis.tilt = 0; 
tomo.rot_axis.corr_area1 = [];
tomo.rot_axis.corr_area2 = [];
tomo.rot_axis.corr_gradient = 0;
tomo.fbp_filter.type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
tomo.fbp_filter.freq_cutoff = 1;
tomo.fbp_filter.padding = 1;
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
% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = '';
write.parfolder = '';
write.subfolder.flatcor = '';
write.subfolder.phase_map = '';
write.subfolder.sino = '';
write.subfolder.reco = '';
write.to_scratch = 0;
write.flatcor = 1;
write.flatcor_shift_cropped = 1;
write.phase_map = 0;
write.sino = 0;
write.phase_sino = 0;
write.reco = 1;
write.float = 1;
write.uint16 = 0;
write.uint8 = 0;
write.post_reconstruction_binning_factor = 2;
write.float_binned = 0;
write.uint16_binned = 0;
write.uint8_binned = 0;
write.flatcor_shift_cropped = 0;
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression.method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression.parameter = [0.02 0.02]; % compression-method specific parameter
% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 0; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 0; % Interactive retrieval to determine regularization parameter
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
poolsize = 0.60; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
% EXPERIMENTAL / NOT YET IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.uint8_segmented = 0; % experimental: threshold segmentation for histograms with 2 distinct peaks: __/\_/\__
compression_method = 'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
compression_parameter = [0.20 0.15]; % compression-method specific parameter
automatic_mode = 0; % Find rotation axis position automatically. NOT IMPLEMENTED!
automatic_mode_coarse = 'entropy'; % NOT IMPLEMENTED!
automatic_mode_fine = 'iso-grad'; % NOT IMPLEMENTED!

% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003440/raw/';

scan_path = [raw_path 'syn01_52L_PEEK_12w']; ADD

scan_path = [raw_path 'syn02_52L_PEEK_12w_ccd']; ADD

scan_path = [raw_path 'syn03_52L_PEEK_12w_dcm_cmos_recotest']; ADD

scan_path = [raw_path 'syn04_52L_PEEK_12w_dcm_cmos_recotest']; ADD

scan_path = [raw_path 'syn05_52L_PEEK_12w_dmm_cmos_recotest']; ADD

scan_path = [raw_path 'syn06_52L_PEEK_12w_dmm_ccd_recotest']; ADD

scan_path = [raw_path 'syn07_cor']; ADD

scan_path = [raw_path 'syn08_cor']; ADD

scan_path = [raw_path 'syn09_cor']; ADD

scan_path = [raw_path 'syn10_cor']; ADD

scan_path = [raw_path 'syn11_cor']; ADD

scan_path = [raw_path 'syn12_cor']; ADD

scan_path = [raw_path 'syn13_cor']; ADD

scan_path = [raw_path 'syn15_cor']; ADD

scan_path = [raw_path 'syn16_39L_PEEK_12w']; ADD

scan_path = [raw_path 'syn17_21R_PEEK_8w']; ADD

scan_path = [raw_path 'syn18_88L_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn19_74R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn20_72R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn21_82R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn22_62L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn23_63L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn24_78L_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn25_69R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn26_69R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn28_69R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn29_78L_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn30_87L_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn31_97R_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn32_99R_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn33_80R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn34_79R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn35_77R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn36_63L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn37_69L_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn38_73R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn39_75L_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn40_69L_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn41_63L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn42_38L_PEEK_8w']; ADD

scan_path = [raw_path 'syn43_38L_PEEK_8w']; ADD

scan_path = [raw_path 'syn44_66L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn45_101BL_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn46_88R_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn47_100AL_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn48_89L_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn49_80L_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn50_99L_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn51_87R_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn52_95R_Mg5Gd_8w']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3-in-1 sample corrosion, CCD
scan_path = [raw_path 'syn60_cor']; raw_bin = 2; rot_axis_offset = -1.5; ADD

% 3-in-1 sample corrosion, CCD
scan_path = [raw_path 'syn61_cor']; raw_bin = 2; rot_axis_offset = -1.5; ADD

% KIT camera
image_correlation.area_width = [0.018 0.038];
im_trafo = 'rot90(im)';
raw_bin = 4; 
bin_before_filtering = 0;
ring_filter.apply = 1;
ring_filter.method = 'wavelet-fft';'jm';
ring_filter.waveletfft.dec_levels = 2:7; 
ring_filter.waveletfft.wname = 'db30'; 'db20';'db25';
ring_filter.waveletfft.sigma = 2.4;
image_correlation.method = 'entropy';
SET_DEFAULT;

% Explant, KIT
scan_path = [raw_path 'syn64_91R_Mg5Gd_4w']; rot_axis_offset = 1.75; ADD

scan_path = [raw_path 'syn65_94R_Mg5Gd_8w']; rot_axis_offset = 1.75; ADD

scan_path = [raw_path 'syn66_cor_Punkt2']; rot_axis_offset = 1.75; ADD

scan_path = [raw_path 'syn67_cor_Punkt1_previous_is_Punkt2']; rot_axis_offset = 2.0; ADD

scan_path = [raw_path 'syn68_cor_Punkt3']; rot_axis_offset = 1.83; ADD

scan_path = [raw_path 'syn69_cor_P1_I'];rot_axis_offset = 1.83;  ADD

scan_path = [raw_path 'syn70_cor_P1_2']; rot_axis_offset = 1.83;  ADD

scan_path = [raw_path 'syn71_cor_P1_3']; rot_axis_offset = 2.0; ADD

scan_path = [raw_path 'syn72_cor_P1_4']; rot_axis_offset = 1.83; ADD

scan_path = [raw_path 'syn73_cor_P1_5']; rot_axis_offset = 1.83; ADD

scan_path = [raw_path 'syn74_cor_P1_6']; rot_axis_offset = 1.5; ADD

scan_path = [raw_path 'syn75_cor_P1_7']; rot_axis_offset = 1.5; ADD

scan_path = [raw_path 'syn76_cor_P2_1']; rot_axis_offset = 1.8; ADD

scan_path = [raw_path 'syn77_cor_P2_2']; rot_axis_offset = 0.8; ADD

scan_path = [raw_path 'syn78_cor_P2_3']; raw_roi = [51 -50]; rot_axis_offset = 1.5; ADD

scan_path = [raw_path 'syn79_cor_P2_4']; raw_roi = []; rot_axis_offset = 1.5 ; ADD

scan_path = [raw_path 'syn80_cor_P2_5']; rot_axis_offset = 1.5 ; ADD

scan_path = [raw_path 'syn81_cor_P2_6']; rot_axis_offset = 1.8 ; ADD

scan_path = [raw_path 'syn82_cor_P2_7']; rot_axis_offset = 1.5 ; ADD

%% Check Rot axis

scan_path = [raw_path 'syn83_cor_P3_1']; ADD; 

scan_path = [raw_path 'syn84_cor_P3_2']; ADD

scan_path = [raw_path 'syn85_cor_P3_3']; ADD

scan_path = [raw_path 'syn86_cor_P3_4']; ADD

scan_path = [raw_path 'syn87_cor_P3_5']; ADD

scan_path = [raw_path 'syn88_cor_P3_6']; ADD

scan_path = [raw_path 'syn89_cor_P3_7']; ADD

scan_path = [raw_path 'syn90_70L']; ADD

scan_path = [raw_path 'syn91_72L']; ADD

scan_path = [raw_path 'syn92_93R']; ADD

scan_path = [raw_path 'syn93_76L']; ADD

scan_path = [raw_path 'syn94_96R']; ADD

scan_path = [raw_path 'syn95_69L']; ADD

scan_path = [raw_path 'syn96_82L']; ADD

scan_path = [raw_path 'syn97_74L']; ADD

scan_path = [raw_path 'syn98_77L']; ADD

scan_path = [raw_path 'syn99_43R']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
