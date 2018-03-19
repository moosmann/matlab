function p05_reco_loop_script_synchroload2017nov2_000( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
% Created on 12-Dec-2017 by moosmanj

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

% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = []; % [y0 y1] vertical roi.  skips first raw_roi(1)-1 lines, reads until raw_roi(2). Not working for *.raw data where images are flipped.
raw_bin = 2; % projection binning factor: 1, 2, or 4
bin_before_filtering = 1; % Binning before pixel filtering is applied, much faster but less effective filtering
proj_range = 1; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = 1; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
dark_FiltPixThresh = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
ref_FiltPixThresh = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
proj_FiltPixThresh = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
correlation_method = 'none';'ssim-ml';'entropy';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';'none';
% 'ssim-ml' : Matlab's structural similarity index (SSIM), includes Gaussian smoothing
% 'ssim' : own implementation of SSIM, smoothing not yet implemented
% 'entropy' : entropy measure of proj over flat
% 'cov' : cross covariance
% 'corr' : cross correlation = normalized cross covariance
% 'std' : standard deviation of proj over flat
% 'diff': difference of proj and flat
% 'shift': computes relative shift from peak of cross-correlation map
% 'none' : no correlation, use median flat
corr_shift_max_pixelshift = 0.25; % maximum pixelshift allowed for 'shift'-correlation method: if 0 use the best match (i.e. the one with the least shift), if > 0 uses all flats with shifts smaller than corr_shift_max_pixelshift
corr_num_flats = 1; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
ring_current_normalization = 1; % normalize flat fields and projections by ring current
flat_corr_area1 = [1 floor(100/raw_bin)];%[0.98 1];% correlation area: index vector or relative/absolute position of [first pix, last pix]
flat_corr_area2 = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
decimal_round_precision = 2; % precision when rounding pixel shifts
ring_filter.apply = 1; % ring artifact filter (only for scans without wiggle di wiggle)
ring_filter.apply_before_stitching = 1; % ! Consider when phase retrieval is applied !
ring_filter.method = 'wavelet-fft';'jm';
ring_filter.waveletfft.dec_levels = 2:6; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft.wname = 'db25';'db30'; % wavelet type for 'wavelet-fft'
ring_filter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter.jm.median_width = 11; % [3 11 21 31 39];
% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; %
% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_tomo = 1; % run tomographic reconstruction
vol_shape = []; % shape of reconstruction volume. in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
vol_size = []; % 6-component vector [xmin xmax ymin ymax zmin zmax]. if empty, volume is centerd within vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs.
rot_angle_full = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
rot_angle_offset = pi; % global rotation of reconstructed volume
rot_axis_offset = []; % if empty use automatic computation
rot_axis_pos = []; % if empty use automatic computation. either offset or pos has to be empty. can't use both
rot_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
rot_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
rot_corr_gradient = 0; % use gradient of intensity maps if signal variations are too weak to correlate projections
rot_axis_tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
fbp_filter_type = 'Ram-Lak';
fpb_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
fbp_filter_padding_method = 'symmetric';
butterworth_filter = 0; % use butterworth filter in addition to FBP filter
butterworth_order = 1;
butterworth_cutoff_frequ = 0.9;
astra_pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted
take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '';
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.flatcor = 0; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if do_tomo=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0;
write.uint8 = 0;
reco_bin = 2; % binning factor of reconstructed volume if binned volumes are saved
write.float_binned = 0; % binned single precision (32-bit float) tiff
write.uint16_binned = 0;
write.uint8_binned = 0;
write.uint8_segmented = 0; % experimental: threshold segmentation for histograms with 2 distinct peaks: __/\_/\__
compression_method = 'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
compression_parameter = [0.20 0.15]; % compression-method specific parameter
% dynamic range is compressed s.t. new dynamic range assumes
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = compression_parameter, eg. [-0.01 1]
% 'std' : NUM = compression_parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = compression_parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
parfolder = '';%sprintf( 'cor_%s', correlation_method);''; % parent folder for 'reco', 'sino', 'phase', and 'flat_corrected'
subfolder_flatcor = ''; % subfolder in 'flat_corrected'
subfolder_phase_map = ''; % subfolder in 'phase_map'
subfolder_sino = ''; % subfolder in 'sino'
subfolder_reco = ''; % subfolder in 'reco'
% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 1; % show images and plots during reconstruction
interactive_determination_of_rot_axis = 1; % reconstruct slices with different rotation axis offsets
interactive_determination_of_rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
lamino = 0; % find laminography tilt instead camera rotation
fixed_tilt = 0; % fixed other tilt
slice_number = 0.5; % slice number, default: 0.5. if in [0,1): relative, if in (1, N]: absolute
% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
poolsize = 0.80; % number of workers used in a parallel pool. if > 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used
link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
automatic_mode = 0; % Find rotation axis position automatically
automatic_mode_coarse = 'entropy'; %
automatic_mode_fine = 'iso-grad';

% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/raw/';

scan_path = [raw_path 'syn134_28R_PEEK_8w']; ADD

scan_path = [raw_path 'syn135_27R_PEEK_8w']; ADD

scan_path = [raw_path 'syn136_95L_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn137_96L_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn138_29R_PEEK_8w']; ADD

scan_path = [raw_path 'syn139_48L_PEEK_12w']; ADD

scan_path = [raw_path 'syn140_61L_PEEK_12w']; ADD

scan_path = [raw_path 'syn141_90L_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn142_35L_PEEK_8w']; ADD

scan_path = [raw_path 'syn143_101BR_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn144_26R_PEEK_8w']; ADD

scan_path = [raw_path 'syn145_58L_Mg_12_cmos_test']; ADD

scan_path = [raw_path 'syn146_58L_Mg_12_cmos_test']; ADD

scan_path = [raw_path 'syn148_58L_Mg_12_cmos_test']; ADD

scan_path = [raw_path 'syn149_58L_Mg_12_cmos_test']; ADD

scan_path = [raw_path 'syn150_58L_Mg_12_000']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_000']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_001']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_002']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_003']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_004']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_005']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_006']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_007']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_008']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_009']; ADD

scan_path = [raw_path 'syn151_58L_Mg_12_011']; ADD

scan_path = [raw_path 'syn152_58L_Mg_12_cmos_test']; ADD

scan_path = [raw_path 'syn153_58L_Mg_12_cmos_test']; ADD

scan_path = [raw_path 'syn154_58L_Mg_12_cmos_test']; ADD

% empty aborted
scan_path = [raw_path 'syn155_berit']; ADD

raw_roi = [1100 2401];
raw_bin = 2;
vol_shape = [0.8 0.8 0.25];
vol_size = [-0.8 0.8 -0.8 0.8 -0.25 0.25];
fbp_filter_type = 'linear';
correlation_method = 'none';
scan_path = [raw_path 'syn156_berit']; ADD

raw_roi = [1100 2401];
raw_bin = 2;
vol_shape = [1.6 1.6 0.5];
vol_size = [-0.8 0.8 -0.8 0.8 -0.5 0.5];
fbp_filter_type = 'Ram-Lak';
butterworth_filter = 1;
correlation_method = 'ssim-ml';
scan_path = [raw_path 'syn156_berit']; ADD

scan_path = [raw_path 'syn157_cor_mg5gd410s_mg10gd401s_png1p_23']; ADD

scan_path = [raw_path 'syn158_cor_mg5gd414s_mg10gd405s_png2p_23']; ADD

scan_path = [raw_path 'syn159_cor_mg5gd430s_mg10gd413s_png9p_23']; ADD

scan_path = [raw_path 'syn160_cor_mg5gd432s_mg10gd418s_png10p_23']; ADD

scan_path = [raw_path 'syn161_cor_mg5gd413s_mg10gd409s_mg5gd3p_23']; ADD

scan_path = [raw_path 'syn162_cor_mg5gd416s_mg10gd410s_mg5gd7p_23']; ADD

scan_path = [raw_path 'syn163_cor_mg5gd428s_mg10gd408s_mg5gd8p_23']; ADD

scan_path = [raw_path 'syn164_cor_mg5gd434s_mg5gd1p_23']; ADD

scan_path = [raw_path 'syn165_58L_Mg10Gd_12w_000']; ADD

scan_path = [raw_path 'syn165_58L_Mg10Gd_12w_001']; ADD

scan_path = [raw_path 'syn165_58L_Mg10Gd_12w_002']; ADD

scan_path = [raw_path 'syn165_58L_Mg10Gd_12w_003']; ADD

scan_path = [raw_path 'syn165_58L_Mg10Gd_12w_004']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_000']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_001']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_002']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_003']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_004']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_005']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_006']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_007']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_008']; ADD

scan_path = [raw_path 'syn166_104R_Mg10Gd_4w_009']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
