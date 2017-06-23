function p05_reco_ehh( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
% Caution: all parameters not defined within the file are taken from
% 'p05_reco'
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
% Written by Julian Moosmann, 2017-06-2, last modification: 2017-06-05
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

visualOutput = 0;
interactive_determination_of_rot_axis = 0;
interactive_determination_of_rot_axis_tilt = 0;
raw_roi = [];
scan_path = '';
raw_bin = 2;
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
gpu_ind = [];

% Set default. Allows parameters to be changed before first data set is added
ADD('default')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_path = '/asap3/petra3/gpfs/p05/2017/data/11002839/raw/';
scan_path = [raw_path 'ehh_2017_019_f'];
raw_roi = [201 2400];
rot_axis_offset = 1052.5 * 2 / raw_bin;
proj_range = 1;
ref_range = 1;
excentric_rot_axis = 1;
crop_at_rot_axis = 0;
stitch_projections = 1; 
stitch_method = 'sine';
do_phase_retrieval = 0;
phase_padding = 1; 
visualOutput = 0;
interactive_determination_of_rot_axis = 0;
ADD('r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
