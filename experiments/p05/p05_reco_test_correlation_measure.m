function p05_reco_test_correlation_measure(nums, doreco, print_field)

if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 0;
end
if nargin < 3
    print_field = 'correlation_method';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sets to loop over %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 0;
default.scan_path = '';
default.visualOutput = 0;
default.scan_path = '';
default.raw_bin = 1;
default.excentric_rot_axis = 0;
default.crop_at_rot_axis = 0;
default.stitch_projections = 0; 
default.proj_range = 1; 
default.ref_range = 1; 
default.correlation_method = 'ssim';
default.corr_num_flats = 3;
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
%% 2017 May
raw = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';
%nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn11_53R_Mg5Gd_12w_load_broken'];

nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_00'];
para(nn).correlation_method =  'ssim-ml';
para(nn).parfolder = 'ssim-ml';

nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_00'];
para(nn).correlation_method =  'ssim';
para(nn).parfolder = 'ssim-LfromROI';

nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_00'];
para(nn).correlation_method =  'ssim';
para(nn).parfolder = 'ssim-L1';

nn = nn + 1;para(nn) = default; para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_00'];
para(nn).correlation_method =  'entropy';
para(nn).parfolder = 'entropy';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( nums, doreco, print_field, para)
