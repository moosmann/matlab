function p05_reco_ivory(nums, doreco)
if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 0;
end
if nargin < 3
    print_field = 'parfolder';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sets to loop over %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003135/raw/';
nn = 0;
default.interactive_determination_of_rot_axis = 0;
default.visualOutput = 0;
default.scan_path = '';
default.raw_bin = 2;
default.excentric_rot_axis = 1;
default.crop_at_rot_axis = 0;
default.stitch_projections = 1; 
default.proj_range = 1; 
default.ref_range = 1; 
default.correlation_method =  'diff';
default.do_phase_retrieval = 1;
default.phase_retrieval_method = 'tie';'qp';'qpcut';
default.phase_retrieval_reg_par = 2; 
default.phase_retrieval_bin_filt = 0.15; 
default.phase_retrieval_cutoff_frequ = 1 * pi; 
default.phase_padding = 1; 
default.do_tomo = 1;
default.vol_shape = [];
default.vol_size = [];
default.ring_filter = 1;
default.rot_axis_offset = [];
default.rot_axis_tilt = 0;
default.fbp_filter_padding = 1;
default.parfolder = 'jm';
default.write_to_scratch = 0;
default.write_flatcor = 1;
default.write_phase_map = 0; 
default.write_sino = 0; 
default.write_sino_phase = 0; 
default.write_reco = 1; 
default.write_float = 1; 
default.write_float_binned = 1; 
default.write_8bit = 1;
default.write_8bit_binned = 1;
default.write_16bit = 0; 

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw_path 'ivo_trans1_006'];
para(nn).rot_axis_offset = 875 / para(nn).raw_bin;
para(nn).phase_retrieval_method = 'tie';
para(nn).phase_retrieval_reg_par = 2; 
para(nn).phase_padding = 1; 

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw_path 'ivo_trans1_006'];
para(nn).rot_axis_offset = 875 / para(nn).raw_bin;
para(nn).phase_retrieval_method = 'qpcut';
para(nn).phase_retrieval_reg_par = 2; 
para(nn).phase_retrieval_bin_filt = 0.15; 
para(nn).phase_retrieval_cutoff_frequ = 1 * pi; 
para(nn).fbp_filter_padding = 0;
para(nn).phase_padding = 0; 

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw_path 'ivo_trans1_006'];
para(nn).raw_bin = 1;
para(nn).rot_axis_offset = 875 / para(nn).raw_bin;
para(nn).phase_retrieval_method = 'qpcut';
para(nn).phase_retrieval_reg_par = 2; 
para(nn).phase_retrieval_bin_filt = 0.15; 
para(nn).phase_retrieval_cutoff_frequ = 1 * pi; 
para(nn).fbp_filter_padding = 0;
para(nn).phase_padding = 0; 
para(nn).vol_shape = [2403, 2403, 1528];
vol_shape = para(nn).vol_shape;
para(nn).vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw_path 'ivo_trans1_006'];
para(nn).raw_bin = 1;
para(nn).rot_axis_offset = 875 / para(nn).raw_bin;
para(nn).phase_retrieval_method = 'tie';
para(nn).phase_retrieval_reg_par = 3; 
para(nn).fbp_filter_padding = 0;
para(nn).phase_padding = 0; 
para(nn).vol_shape = [2403, 2403, 1528];
vol_shape = para(nn).vol_shape;
para(nn).vol_size = [-vol_shape(1), vol_shape(1), -vol_shape(2), vol_shape(2), -vol_shape(3), vol_shape(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( nums, doreco, print_field, para)
