function p05_reco_test_correlation_measure(SUBSETS, RUN_RECO, PRINT_PARAMETERS )
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

scan_path = '';
visualOutput = 0;
scan_path = '';
raw_bin = 1;
excentric_rot_axis = 0;
crop_at_rot_axis = 0;
stitch_projections = 0; 
proj_range = 1; 
ref_range = 1; 
correlation_method = 'ssim';
corr_num_flats = 3;
do_phase_retrieval = 0;
phase_retrieval_method = 'tie';'qp';'qpcut';
phase_retrieval_reg_par = 2.5; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
do_tomo = 1;
ring_filter = 1;
rot_axis_offset = [];
rot_axis_tilt = [];
parfolder = '';
write_to_scratch = 0;
write_flatcor = 1;
write_phase_map = 0; 
write_sino = 0; 
write_sino_phase = 0; 
write_reco = 1; 
write_float = 1; 
write_float_binned = 1; 
write_8bit = 0;
write_8bit_binned = 0;
write_16bit = 0; 

ADD_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2017 May
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';
%scan_path = [raw_path 'syn11_53R_Mg5Gd_12w_load_broken'];
ADD_DATA_SET

scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
correlation_method =  'ssim-ml';
parfolder = 'ssim-ml';
ADD_DATA_SET

scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
correlation_method =  'ssim';
parfolder = 'ssim-LfromROI';
ADD_DATA_SET

scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
correlation_method =  'ssim';
parfolder = 'ssim-L1';
ADD_DATA_SET

scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
correlation_method =  'entropy';
parfolder = 'entropy';
ADD_DATA_SET

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS, PARAMETER_CELL)
