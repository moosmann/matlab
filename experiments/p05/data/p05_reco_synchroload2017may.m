function p05_reco_synchroload2017may( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2017 May
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';

scan_path = [raw_path 'syn01_48L_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn01_48L_PEEK_12w_c'];ADD
scan_path = [raw_path 'syn02_46R_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn02_46R_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn03_12R_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn03_12R_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn04_30R_PEEK_8w_a'];ADD
scan_path = [raw_path 'syn04_30R_PEEK_8w_b'];ADD
scan_path = [raw_path 'syn05_41R_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn11_53R_Mg5Gd_12w_load_broken'];ADD

raw_roi = [287 1486];
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_02'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_04'];
ref_range = [1:135, 137:162];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_06'];
ref_range = 1;ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_08'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_10'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_12'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_14'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_16'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_18'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_20'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_22'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_24_noload'];ADD

raw_roi = [1211 2410];
scan_path = [raw_path 'syn14_48L_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn14_48L_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn15_29R_PEEK_8w_a'];ADD
scan_path = [raw_path 'syn15_29R_PEEK_8w_b'];ADD
scan_path = [raw_path 'syn16_43R_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn16_43R_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn17_25R_PEEK_8w_a'];ADD('r')

%% Radiography after load increase before tomography
raw_roi = [];
do_tomo = 0;
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_01'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_03'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_05'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_07'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_09'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_11'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_13'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_15'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_17'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_19'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_21'];ADD('r')

%% Tomography CPD straw 2
raw_roi = [ 211 1420];
scan_path = [raw_path 'syn22_77L_Mg5Gd_8w_a'];ADD
scan_path = [raw_path 'syn22_77L_Mg5Gd_8w_b'];ADD
scan_path = [raw_path 'syn22_80L_Mg5Gd_8w_a'];ADD
scan_path = [raw_path 'syn22_80L_Mg5Gd_8w_b'];ADD
scan_path = [raw_path 'syn22_87L_Mg5Gd_4w_a'];ADD
scan_path = [raw_path 'syn22_87L_Mg5Gd_4w_b'];ADD
scan_path = [raw_path 'syn22_88R_Mg5Gd_4w_a'];ADD
scan_path = [raw_path 'syn22_88R_Mg5Gd_4w_b'];ADD
scan_path = [raw_path 'syn22_99L_Mg5Gd_4w_a'];ADD
scan_path = [raw_path 'syn22_99L_Mg5Gd_4w_b'];ADD('r')

% CPD straw II:top
scan_path = [raw_path 'syn23_28R_PEEK_8w_a'];ADD
scan_path = [raw_path 'syn23_28R_PEEK_8w_b'];ADD('r')


%% TEST SECTION
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
raw_roi = [287 1486];
raw_bin = 1;
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';
phase_retrieval_reg_par = 2.5; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
write_8bit_binned = 1;
ADD

phase_retrieval_reg_par = 1.5; ADD
phase_retrieval_reg_par = 3.5; ADD
phase_retrieval_method = 'qp';phase_retrieval_reg_par = 2.5; ADD
phase_retrieval_method = 'qpcut';ADD('r')

scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
raw_roi = [287 1486];
write_float = 1; 
write_float_binned = 1; 
write_16bit = 1;
write_16bit_binned = 1;
write_8bit = 1;
write_8bit_binned = 1;
compression_method = 'histo';
compression_parameter = [0.05 0.05];
parfolder = 'test';
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
