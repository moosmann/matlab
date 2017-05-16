function p05_reco_synchroload2017may(nums, doreco)

if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 1;
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
default.write_8bit = 1;
default.write_8bit_binned = 1;
default.write_16bit = 0; 

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2017 May
raw = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';

nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn01_48L_PEEK_12w_b'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn01_48L_PEEK_12w_c'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn02_46R_PEEK_12w_a'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn02_46R_PEEK_12w_b'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn03_12R_PEEK_12w_a'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn03_12R_PEEK_12w_b'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn04_30R_PEEK_8w_a'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn04_30R_PEEK_8w_b'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn05_41R_PEEK_12w_a'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn11_53R_Mg5Gd_12w_load_broken'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_00'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_02'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_04'];para(nn).ref_range = [1:135, 137:162];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_06'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_08'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_10'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_12'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_14'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_16'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_18'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_20'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_22'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_24_noload'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn14_48L_PEEK_12w_a'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn14_48L_PEEK_12w_b'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn15_29R_PEEK_8w_a'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn15_29R_PEEK_8w_b'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn16_43R_PEEK_12w_a'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn16_43R_PEEK_12w_b'];
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn17_25R_PEEK_8w_a'];

% Radio after load increase before tomography
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_01'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_03'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_05'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_07'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_09'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_11'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_13'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_15'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_17'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_19'];para(nn).do_tomo = 0;
nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_21'];para(nn).do_tomo = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( nums, doreco, para)
