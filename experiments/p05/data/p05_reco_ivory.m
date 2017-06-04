function p05_reco_ivory(SUBSETS, RUN_RECO, PRINT_PARAMETERS )
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

raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003135/raw/';
interactive_determination_of_rot_axis = 0;
visualOutput = 0;
scan_path = '';
raw_bin = 2;
excentric_rot_axis = 1;
crop_at_rot_axis = 0;
stitch_projections = 1; 
proj_range = 1; 
ref_range = 1; 
correlation_method =  'diff';
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';'qp';'qpcut';
phase_retrieval_reg_par = 2; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
do_tomo = 1;
vol_shape = [];
vol_size = [];
ring_filter = 1;
rot_axis_offset = [];
rot_axis_tilt = 0;
fbp_filter_padding = 1;
parfolder = 'jm';
write_to_scratch = 0;
write_flatcor = 1;
write_phase_map = 0; 
write_sino = 0; 
write_sino_phase = 0; 
write_reco = 1; 
write_float = 1; 
write_float_binned = 1; 
write_8bit = 1;
write_8bit_binned = 1;
write_16bit = 0; 

ADD('d')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scan_path = [raw_path 'ivo_trans1_006'];
rot_axis_offset = 875 / raw_bin;
phase_retrieval_method = 'tie';
phase_retrieval_reg_par = 2; 
phase_padding = 1; 
ADD

scan_path = [raw_path 'ivo_trans1_006'];
rot_axis_offset = 875 / raw_bin;
phase_retrieval_method = 'qpcut';
phase_retrieval_reg_par = 2; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
fbp_filter_padding = 0;
phase_padding = 0; 
ADD


scan_path = [raw_path 'ivo_trans1_006'];
raw_bin = 1;
rot_axis_offset = 875 / raw_bin;
phase_retrieval_method = 'qpcut';
phase_retrieval_reg_par = 2; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
fbp_filter_padding = 0;
phase_padding = 0; 
vol_shape = [2403, 2403, 1528];
vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];
ADD


scan_path = [raw_path 'ivo_trans1_006'];
raw_bin = 1;
rot_axis_offset = 875 / raw_bin;
phase_retrieval_method = 'tie';
phase_retrieval_reg_par = 3; 
fbp_filter_padding = 0;
phase_padding = 0; 
vol_shape = [2403, 2403, 1528];
vol_size = [-vol_shape(1), vol_shape(1), -vol_shape(2), vol_shape(2), -vol_shape(3), vol_shape(3)];
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
