function p05_reco_loop_hnee_tension_wood_000( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
% Created on 02-Jul-2018 by moosmanj

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

read_flatcor = 0;
energy = []; 
sample_detector_distance = []; 
eff_pixel_size = []; 
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = [301 3300];
%raw_roi = -1; 
raw_bin = 2; 
bin_before_filtering(1) = 0; 
excentric_rot_axis = 0;
crop_at_rot_axis = 0; 
stitch_projections = 0; 
proj_range = []; 
ref_range = []; 
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixe
ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';
image_correlation.num_flats = 1; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
image_correlation.area_width = [0 0.02];%[0.98 1];% correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.area_height = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.shift.max_pixelshift = 0.25;
ring_filter.apply = 0; 
ring_filter.apply_before_stitching = 0; 
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft.dec_levels = 2:5; 
ring_filter.waveletfft.wname = 'db25';'db30';
ring_filter.waveletfft.sigma = 2.4;
ring_filter.jm.median_width = 11; 
strong_abs_thresh = 1; 
decimal_round_precision = 2; 
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 1; 
phase_retrieval.apply_before = 1; 
phase_retrieval.post_binning_factor = 1;
phase_retrieval.method = 'tie';'qpcut'; 
phase_retrieval.reg_par = 1.5; 
phase_retrieval.bin_filt = 0.15;
phase_retrieval.cutoff_frequ = 2 * pi;
phase_retrieval.padding = 1;
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; 
tomo.run_interactive_mode = 0;
tomo.reco_mode = '3D'; 
tomo.vol_size = [-0.5 0.5 -0.5 0.5 -0.5 0.5];
tomo.vol_shape = [];
tomo.rot_axis.offset = [];
tomo.rot_axis.position = []; 
tomo.rot_angle.full_range = []; 
tomo.rot_angle.offset = pi / 2; 
tomo.rot_axis.tilt = 0; 
tomo.rot_axis.corr_area1 = [];
tomo.rot_axis.corr_area2 = [];
tomo.rot_axis.corr_gradient = 0;
tomo.fbp_filter.type = 'linear';
tomo.fbp_filter.freq_cutoff = 1;
tomo.fbp_filter.padding = 1; 
tomo.fbp_filter.padding_method = 'symmetric';
tomo.butterworth_filter.apply = 1; % use butterworth filter in addition to FBP filter
tomo.butterworth_filter.order = 1;
tomo.butterworth_filter.frequ_cutoff = 0.95;
tomo.astra_pixel_size = 1; 
tomo.take_neg_log = []; 
tomo.algorithm = 'fbp';'sirt'; 'cgls';
tomo.iterations = 40; % for 'sirt' or 'cgls'.
tomo.sirt.MinConstraint = []; 
tomo.sirt.MaxConstraint = []; 
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = '';
write.to_scratch = 0; 
write.parfolder = '';
write.subfolder.flatcor = '';
write.subfolder.phase_map = '';
write.subfolder.sino = ''; 
write.subfolder.reco = ''; 
write.flatcor = 1; 
write.phase_map = 0;
write.sino = 0; 
write.phase_sino = 0; 
write.reco = 1; 
write.float = 1; 
write.uint16 = 0;
write.uint8 = 0; 
% Optionally save binned reconstructions
write.float_binned = 1;
write.uint16_binned = 0;
write.uint8_binned = 0; 
write.reco_binning_factor = 2;
write.compression.method = 'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression.parameter = [0.20 0.15]; % compression-method specific parameter
write.uint8_segmented = 0;
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 0; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 0; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
poolsize = 0.60; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all

% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2018/data/11004450/raw/';

%image_correlation.method = 'none';
%raw_roi = [2301 2600];
%write.reco = 0;
%tomo.run = 0;

%scan_path = [raw_path 'hnee01_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee02_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee03_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee04_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee05_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee06_kie_sh_ts_000N']; ADD % some mismatch
%raw_roi = [301 3300];
tomo.rot_axis.offset = 2 * -784 / raw_bin;
scan_path = [raw_path 'hnee07_kie_sh_ts_000N']; ADD

%raw_roi = [301 3300];
tomo.rot_axis.offset = 2 * -783.5 / raw_bin;
scan_path = [raw_path 'hnee08_kie_sh_ts_000']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_001']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_002']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_003']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_004']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_005']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_006']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_007']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_008']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_009']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_010']; ADD

scan_path = [raw_path 'hnee10_kie_fh_ts_0p5N']; ADD
scan_path = [raw_path 'hnee11_kie_fh_ts_0p5N']; ADD
scan_path = [raw_path 'hnee12_kie_fh_ts_0p5N']; ADD
scan_path = [raw_path 'hnee13_kie_fh_ts_0p5N']; ADD
scan_path = [raw_path 'hnee14_kie_fh_ts_0p5N']; ADD

% not working
% scan_path = [raw_path 'hnee15_kie_fh_ts_000']; ADD
%scan_path = [raw_path 'hnee16_kie_fh_ts_000']; ADD
%scan_path = [raw_path 'hnee17_kie_fh_ts_000']; ADD

tomo.rot_axis.offset = 2 * -10.0 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_000']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_001']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_002']; ADD
tomo.rot_axis.offset = 2 * -10.85 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_003']; ADD
tomo.rot_axis.offset = 2 * -11.0 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_004']; ADD
scan_path = [raw_path 'hnee18_pappel_tensionWood_005']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_006']; ADD
tomo.rot_axis.offset = 2 * -11.0 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_007']; ADD
tomo.rot_axis.offset = 2 * -11.25 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_008']; ADD

tomo.rot_axis.offset = 2 * -10.85 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_000']; ADD
tomo.rot_axis.offset = 2 * -10.3 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_001']; ADD
tomo.rot_axis.offset = 2 * -10.35 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_002']; ADD
tomo.rot_axis.offset = 2 * -10.35 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_003']; ADD 
% bad reco
scan_path = [raw_path 'hnee19_pappel_oppositeWood_004']; ADD
tomo.rot_axis.offset = 2 * -10.6 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_005']; ADD
tomo.rot_axis.offset = 2 * -10.7 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_006']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_007']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_008']; ADD

% Movement
tomo.rot_axis.offset = 2 * 8.5 / raw_bin;
scan_path = [raw_path 'hnee20_pappel_tensionWood_000']; ADD
tomo.rot_axis.offset = 2 * 8.5 / raw_bin;
scan_path = [raw_path 'hnee20_pappel_tensionWood_001']; ADD
tomo.rot_axis.offset = 2 * 8.25 / raw_bin;
scan_path = [raw_path 'hnee20_pappel_tensionWood_002']; ADD
tomo.rot_axis.offset = 2 * 8.0 / raw_bin;
scan_path = [raw_path 'hnee20_pappel_tensionWood_003']; ADD

tomo.rot_axis.offset = 2 * 7.5 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_000']; ADD
tomo.rot_axis.offset = 2 * 9.0 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_001']; ADD
tomo.rot_axis.offset = 2 * 8.85 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_002']; ADD
tomo.rot_axis.offset = 2 * 8.9 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_003']; ADD
tomo.rot_axis.offset = 2 * 8.9 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_004']; ADD
tomo.rot_axis.offset = 2 * 8.8 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_005']; ADD
tomo.rot_axis.offset = 2 * 9.05 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_006']; ADD

tomo.rot_axis.offset = 2 * 9.1 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_000']; ADD
tomo.rot_axis.offset = 2 * 9.1 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_001']; ADD
tomo.rot_axis.offset = 2 * 9.1 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_002']; ADD
tomo.rot_axis.offset = 2 * 10.3 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_003']; ADD
tomo.rot_axis.offset = 2 * 10.3 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_004']; ADD
tomo.rot_axis.offset = 2 * 10.39 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_005']; ADD

tomo.rot_axis.offset = 2 * 9.5 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_000']; ADD
tomo.rot_axis.offset = 2 * 9.9 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_001']; ADD
tomo.rot_axis.offset = 2 * 9.95 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_002']; ADD
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_003']; ADD
tomo.rot_axis.offset = 2 * 10.05 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_004']; ADD
% No top up beam current drops from 92 to 82
tomo.rot_axis.offset = 2 * 10.10 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_005']; ADD
% No top up beam current drops from 72? to 
tomo.rot_axis.offset = 2 * 10.05 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_006']; ADD
% No top up beam current drops from 66 to 60
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_007']; ADD
% No top up beam current drops from 58 to 55
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_008']; ADD
% No top up beam current drops from 52 to 48
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_009']; ADD
% No top up beam current drops from 47 to 44
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_010']; ADD
% No top up beam current drops from 42 to 40
tomo.rot_axis.offset = 2 * 9.9 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_011']; ADD
% No top up beam current drops from 38 to 36
tomo.rot_axis.offset = 2 * 9.95 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_012']; ADD
% No top up beam current drops from 35 to 33
%% CHECK
tomo.rot_axis.offset = 2 * 10 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_013']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
