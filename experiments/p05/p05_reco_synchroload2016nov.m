function p05_reco_synchroload2016nov(nums, doreco, print_field)

if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 1;
end
if nargin < 3
    print_field = '';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sets to loop over %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 0;
default.visualOutput = 0;
default.scan_path = '';
default.raw_bin = 1;
default.excentric_rot_axis = 0;
default.crop_at_rot_axis = 0;
default.stitch_projections = 0; 
default.proj_range = 1; 
default.ref_range = 1; 
default.correlation_method =  'ssim-ml';
default.do_phase_retrieval = 0;
default.phase_retrieval_method = 'tie';'qp';'qpcut';
default.phase_retrieval_reg_par = 2.5; 
default.phase_retrieval_bin_filt = 0.15; 
default.phase_retrieval_cutoff_frequ = 1 * pi; 
default.phase_padding = 1; 
default.do_tomo = 1;
default.ring_filter = 1;
default.rot_axis_offset = [];
default.rot_axis_tilt = [];
default.parfolder = '';
default.write_to_scratch = 0;
default.write_flatcor = 1;
default.write_phase_map = 0; 
default.write_sino = 0; 
default.write_sino_phase = 0; 
default.write_reco = 1; 
default.write_float = 1; 
default.write_float_binned = 1; 
default.write_8bit = 0;
default.write_8bit_binned = 0;
default.write_16bit = 0; 

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/';

% corroded screw
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'mah_01'];
para(nn).rot_axis_offset = -135.75 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corroded screw
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'mah_02'];
para(nn).rot_axis_offset = -135.75 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corroded screw
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'mah_03'];
para(nn).rot_axis_offset = -135.75 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% implant fresh
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'mah_04'];
para(nn).excentric_rot_axis = 1;
para(nn).crop_at_rot_axis = 1;
para(nn).rot_axis_offset = 628 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;
%para(nn).vol_shape = [2155 2155 1050];

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_05'];
para(nn).rot_axis_offset = 2 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_06_Mg10G_004'];
para(nn).rot_axis_offset = 5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% implant fresh formalin
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_07_bone_in_formalin'];
para(nn).rot_axis_offset = 2 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corrosion cell
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_08_corrosion_cell_A'];
para(nn).rot_axis_offset = -3.5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corrosion cell
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_08_corrosion_cell_B'];
para(nn).rot_axis_offset = -2 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_10_13R_top'];
para(nn).excentric_rot_axis = 1;
para(nn).rot_axis_offset = 1078 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.00267; % about -.15 degrees

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw '/mah_10_13R_bottom'];
para(nn).rot_axis_offset = 539;
para(nn).rot_axis_tilt = -0.0023;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_11_20R_top'];
para(nn).rot_axis_offset = 538.5;
para(nn).rot_axis_tilt = -0.0024;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_11_20R_bottom'];
para(nn).rot_axis_offset = 538.5;
para(nn).rot_axis_tilt = -0.0024;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_15_57R'];
para(nn).excentric_rot_axis = 0;
para(nn).crop_at_rot_axis = 0;
para(nn).rot_axis_offset = -2.5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.0028;
para(nn).do_phase_retrieval = 0;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_15_57R'];
para(nn).excentric_rot_axis = 0;
para(nn).crop_at_rot_axis = 0;
para(nn).rot_axis_offset = -2.5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.0028;
para(nn).do_phase_retrieval = 1;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_16_57R_load'];
para(nn).excentric_rot_axis = 0;
para(nn).crop_at_rot_axis = 0;
para(nn).rot_axis_offset = -2.5;
para(nn).rot_axis_tilt = -0.0028;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_17_57R_load_middle'];
para(nn).rot_axis_offset = 88.25;
para(nn).rot_axis_tilt = -0.0025;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_18_57R_load_top'];
para(nn).rot_axis_offset = 88.25;
para(nn).rot_axis_tilt = -0.003;

% 3.4.17
nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_20_4L_bottom'];
para(nn).rot_axis_offset = 12.5;
para(nn).rot_axis_tilt = -0.0029;
para(nn).do_phase_retrieval = 0;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).do_phase_retrieval = 1;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_22_50L_top'];
para(nn).rot_axis_offset = 88.5;
para(nn).rot_axis_tilt = -0.0028;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_23_50L_top'];
para(nn).rot_axis_offset = 5;
para(nn).rot_axis_tilt = -0.004160;
para(nn).do_phase_retrieval = 0;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_24_50L_top_load'];
para(nn).rot_axis_offset = 4.5;
para(nn).rot_axis_tilt = -0.004;
para(nn).do_phase_retrieval = 1;


nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_28_15R_top'];
para(nn).rot_axis_offset = 8.25;
para(nn).rot_axis_tilt = -0.003;
para(nn).do_phase_retrieval = 0;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_30_15R_top_occd125_withoutpaper'];
para(nn).rot_axis_offset = -3.5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;
para(nn).do_phase_retrieval = 1;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_32_15R_top_occd800_withoutpaper'];
para(nn).rot_axis_offset = -79.0 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.002;
para(nn).do_phase_retrieval = 1;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_33_50L_occd400_bottom'];
para(nn).rot_axis_offset = -40 / para(nn).raw_bin;
para(nn).rot_axis_tilt = 0.00137;
para(nn).do_phase_retrieval = 1;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_33_50L_occd400_top'];
para(nn).rot_axis_offset = 5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = 0.00135;
para(nn).do_phase_retrieval = 1;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_35_1R_bottom'];
para(nn).rot_axis_offset = 2 / para(nn).raw_bin;
para(nn).rot_axis_tilt = 0.00158;
para(nn).do_phase_retrieval = 0;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_36_1R_top'];
para(nn).rot_axis_offset = 2 / para(nn).raw_bin;
para(nn).rot_axis_tilt = 0.00158;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_37_10R_bottom'];
para(nn).rot_axis_offset = 0 / para(nn).raw_bin;
para(nn).rot_axis_tilt = 0.00158;

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_38_10R_top'];

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_39_3L_bottom'];

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_40_3L_top'];

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_41_9R_bottom'];

nn = nn + 1;
para(nn) = para(nn - 1);
para(nn).scan_path = [ raw 'mah_42_9R_top'];

% Straw: no proper reco possible due to movment

% corroded screw: movement
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_01'];
para(nn).rot_axis_offset = 2 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corroded screw: movement
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_02'];
para(nn).rot_axis_offset = -1.5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_03'];
para(nn).rot_axis_offset = -1.25 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.0027;
% time-varying bright spots: for nn=1:40,imsc(flat(0+(1:400),0+(1:200),nn)',[000 9000]),pause(1),end

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_04'];
para(nn).rot_axis_offset = -1.75 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_05'];
para(nn).rot_axis_offset = -2.5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.0025;

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_06'];
para(nn).rot_axis_offset = -0 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_2_00'];
para(nn).rot_axis_offset = -0.4 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_2_00'];
para(nn).rot_axis_offset = -0.4 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;
para(nn).raw_bin = 2;
para(nn).parfolder = 'bin2';
para(nn).write_float_binned = 0;
para(nn).write_8bit_binned = 0;

% corroded screw
nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'mah_straw_2_01'];
para(nn).rot_axis_offset = -0 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( nums, doreco, print_field, para)
