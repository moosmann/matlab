function p05_reco_loop_template( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Template function which loops over the data sets given in the 'PARAMETER /
% DATA SETS' section below. The 'DEFAULT PARAMETERS' section defines the
% default paramters. Data / parameter sets are added to the loop using the
% command 'ADD_DATA_SET'. 
%
% Caution: if parameters are changed, they remain changed after the data
% set is added unless 'ADD_DATA_SET(1)' is used which restores all
% parameters given in the 'DEFAULT PARAMETER' section.
%
% Caution: If parameters not set in the DEFAULT PARAMETER section, then
% the values in main reconstruction routine 'p05_reco' are used.
% 
% Copy this dummy function under a new name and add parameter or data sets
% to loop over.
% 
% ARGUMENTS
% SUBSETS : 1D array of integers. subset of data sets to be looped over
% RUN_RECO : bool. default: 0. 0: loops over the subsets but does not start
% reconstructions, 1: start the reconstruction loop.
% PRINT_PARAMETERS : string or cell of strings. parameter to be printed at each
% loop step. useful in combination with RUN_RECO = 0 to check parameter
% setting for the sets to loop over
%
% Written by Julian Moosmann, 2017-06-2, last modification: 2017-06-03
%
% p05_reco_loop_template( SUBSETS, RUN_RECO, PRINT_PARAMETERS)

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

% Example: modify, replace, delete
visualOutput = 0;
interactive_determination_of_rot_axis = 0;
interactive_determination_of_rot_axis_tilt = 0;
raw_roi = [];
scan_path = '';
raw_bin = 1;
excentric_rot_axis = 0;
crop_at_rot_axis = 0;
stitch_projections = 0; 
proj_range = 1; 
ref_range = 1; 
correlation_method =  'ssim-ml';
do_phase_retrieval = 0;
phase_retrieval_method = 'tie';'qp';'qpcut';
phase_retrieval_reg_par = 2.5; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
do_tomo = 1;
fbp_filter_padding = 1;
fbp_filter_freq_cutoff = 1;
ring_filter = 1;
dec_levels = 7;
wname = 'db25';
sigma = 2.4;
rot_axis_offset = [];
rot_axis_tilt = 0.001;
parfolder = '';
write_to_scratch = 0;
write_flatcor = 0;
write_phase_map = 0; 
write_sino = 0; 
write_sino_phase = 0; 
write_reco = 1; 
write_float = 1; 
write_float_binned = 1; 
write_8bit = 0;
write_8bit_binned = 1;
write_16bit = 0; 
subfolder_reco = '';
gpu_ind = 1;

% Set default. Sllows parameters to be changed before first data set is added
ADD_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example: modify, replace, delete

raw = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';

% Define scan path and add data set
scan_path = [raw 'syn01_48L_PEEK_12w_b'];
ADD_DATA_SET;

% Add another data set
scan_path = [raw 'syn01_48L_PEEK_12w_c'];
ADD_DATA_SET();

% Change parameter, add data set, and restore default parameters
scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_04'];
ref_range = [1:135, 137:162];
ADD_DATA_SET(1);

% Changer paramter for subsequent data sets and add data sets
raw_roi = [1211 2410];
scan_path = [raw 'syn14_48L_PEEK_12w_a'];
ADD_DATA_SET();

scan_path = [raw 'syn14_48L_PEEK_12w_b'];
ADD_DATA_SET();

% Add data set and restore defaults
scan_path = [raw 'syn17_25R_PEEK_8w_a'];ADD_DATA_SET(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS, PARAMETER_CELL)
