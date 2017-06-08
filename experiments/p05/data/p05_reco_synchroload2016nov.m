function p05_reco_synchroload2016nov(SUBSETS, RUN_RECO, PRINT_PARAMETERS )
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

visualOutput = 1;
interactive_determination_of_rot_axis = 1;
scan_path = '';
raw_roi = [];
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
ring_filter = 1;
ring_filter_method =  'jm';'wavelet-fft';
ring_filter_median_witdth = 11;
rot_axis_offset = [];
rot_axis_tilt = [];
parfolder = '';
write_to_scratch = 0;
write_flatcor = 0;
write_phase_map = 0; 
write_sino = 0; 
write_sino_phase = 0; 
write_reco = 1; 
write_float = 0; 
write_float_binned = 1; 
write_8bit = 0;
write_8bit_binned = 1;
write_16bit = 0; 
dec_levels = 6;
wname = 'db30';
sigma = 2.4;
gpu_index = [];

ADD_DEFAULT

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/';

% corroded screw
scan_path = [raw_path 'mah_01'];
rot_axis_offset = -135.75 / raw_bin;
rot_axis_tilt = -0.003;
%write_sino = 1;
ADD

% corroded screw
scan_path = [raw_path 'mah_02'];
rot_axis_offset = -135.75 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw
scan_path = [raw_path 'mah_03'];
rot_axis_offset = -135.75 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% implant fresh
scan_path = [raw_path 'mah_04'];
excentric_rot_axis = 1;
crop_at_rot_axis = 1;
rot_axis_offset = 628 / raw_bin;
rot_axis_tilt = -0.003;
%vol_shape = [2155 2155 1050];
ADD('r')

% corroded screw
scan_path = [raw_path 'mah_05'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw
scan_path = [raw_path 'mah_06_Mg10G_004'];
rot_axis_offset = 5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% implant fresh formalin
scan_path = [raw_path 'mah_07_bone_in_formalin'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corrosion cell
scan_path = [raw_path 'mah_08_corrosion_cell_A'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corrosion cell
scan_path = [raw_path 'mah_08_corrosion_cell_B'];
rot_axis_offset = -2 / raw_bin;
rot_axis_tilt = -0.003;
ADD

%
scan_path = [raw_path 'mah_10_13R_top'];
excentric_rot_axis = 1;
rot_axis_offset = 1078 / raw_bin;
rot_axis_tilt = -0.00267; % about -.15 degrees
ADD


scan_path = [ raw_path '/mah_10_13R_bottom'];
rot_axis_offset = 539;
rot_axis_tilt = -0.0023;
ADD

scan_path = [ raw_path 'mah_11_20R_top'];
rot_axis_offset = 538.5;
rot_axis_tilt = -0.0024;
ADD

scan_path = [ raw_path 'mah_11_20R_bottom'];
rot_axis_offset = 538.5;
rot_axis_tilt = -0.0024;
ADD('r')

scan_path = [ raw_path 'mah_15_57R'];
rot_axis_offset = -2.5 / raw_bin;
rot_axis_tilt = -0.0028;
ADD

scan_path = [ raw_path 'mah_15_57R'];
rot_axis_offset = -2.5 / raw_bin;
rot_axis_tilt = -0.0028;
do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_16_57R_load'];
rot_axis_offset = -2.5;
rot_axis_tilt = -0.0028;
ADD

scan_path = [ raw_path 'mah_17_57R_load_middle'];
rot_axis_offset = 88.25;
rot_axis_tilt = -0.0025;
ADD

scan_path = [ raw_path 'mah_18_57R_load_top'];
rot_axis_offset = 88.25;
rot_axis_tilt = -0.003;
ADD

% 3.4.17
scan_path = [ raw_path 'mah_20_4L_bottom'];
rot_axis_offset = 12.5;
rot_axis_tilt = -0.0029;
ADD

do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_22_50L_top'];
rot_axis_offset = 88.5;
rot_axis_tilt = -0.0028;
ADD

scan_path = [ raw_path 'mah_23_50L_top'];
rot_axis_offset = 5;
rot_axis_tilt = -0.004160;
ADD

scan_path = [ raw_path 'mah_24_50L_top_load'];
rot_axis_offset = 4.5;
rot_axis_tilt = -0.004;
do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_28_15R_top'];
rot_axis_offset = 8.25;
rot_axis_tilt = -0.003;
ADD

scan_path = [ raw_path 'mah_29_15R_top_occd125_withpaper'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

scan_path = [ raw_path 'mah_29_15R_top_occd125_withpaper'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_30_15R_top_occd125_withoutpaper'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

scan_path = [ raw_path 'mah_30_15R_top_occd125_withoutpaper'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_32_15R_top_occd800_withoutpaper'];
rot_axis_offset = -79.0 / raw_bin;
rot_axis_tilt = -0.002;
do_phase_retrieval = 1;
ADD

scan_path = [ raw_path 'mah_33_50L_occd400_bottom'];
rot_axis_offset = -40 / raw_bin;
rot_axis_tilt = 0.00137;
do_phase_retrieval = 1;
ADD

scan_path = [ raw_path 'mah_33_50L_occd400_top'];
rot_axis_offset = 5 / raw_bin;
rot_axis_tilt = 0.00135;
do_phase_retrieval = 1;
ADD

scan_path = [ raw_path 'mah_35_1R_bottom'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = 0.00158;
do_phase_retrieval = 0;
ADD

scan_path = [ raw_path 'mah_36_1R_top'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = 0.00158;
ADD

scan_path = [ raw_path 'mah_37_10R_bottom'];
raw_roi = [1 1100];
rot_axis_offset = 0.6 / raw_bin;
rot_axis_tilt = 0.0015;
ADD

scan_path = [ raw_path 'mah_38_10R_top'];
rot_axis_offset = 0.6 / raw_bin;
rot_axis_tilt = 0.0015;
ADD

scan_path = [ raw_path 'mah_39_3L_bottom'];
ADD

scan_path = [ raw_path 'mah_40_3L_top'];
ADD

scan_path = [ raw_path 'mah_41_9R_bottom'];
ADD

scan_path = [ raw_path 'mah_42_9R_top'];
ADD

% Straw: no proper reco possible due to movment

% corroded screw: movement
scan_path = [raw_path 'mah_straw_01'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw: movement
scan_path = [raw_path 'mah_straw_02'];
rot_axis_offset = -1.5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw
scan_path = [raw_path 'mah_straw_03'];
rot_axis_offset = -1.25 / raw_bin;
rot_axis_tilt = -0.0027;
% time-varying bright spots: for nn=1:40,imsc(flat(0+(1:400),0+(1:200),nn)',[000 9000]),pause(1),end

% corroded screw
scan_path = [raw_path 'mah_straw_04'];
rot_axis_offset = -1.75 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw
scan_path = [raw_path 'mah_straw_05'];
rot_axis_offset = -2.5 / raw_bin;
rot_axis_tilt = -0.0025;

% corroded screw
scan_path = [raw_path 'mah_straw_06'];
rot_axis_offset = -0 / raw_bin;
rot_axis_tilt = -0.003;

% corroded screw
scan_path = [raw_path 'mah_straw_2_00'];
rot_axis_offset = -0.4 / raw_bin;
rot_axis_tilt = -0.003;
write_sino = 1; 

% corroded screw
scan_path = [raw_path 'mah_straw_2_01'];
rot_axis_offset = -0 / raw_bin;
rot_axis_tilt = -0.003;
write_sino = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
