function p05_reco_synchroload(nums, doreco)

if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 0;
end
if nargin < 3
    print_field = '';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sets to loop over %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 0;
default.raw_roi = [];
default.scan_path = '';
default.raw_bin = 1;
default.excentric_rot_axis = 0;
default.crop_at_rot_axis = 0;
default.stitch_projections = 0; 
default.proj_range = 1; 
default.ref_range = 1; 
default.correlation_method =  'ssim-ml';
default.do_phase_retrieval = 0;
default.do_tomo = 1;
default.fbp_filter_padding = 1;
default.ring_filter = 1;
default.ring_filter_method = 'jm';
default.ring_filter_median_width = 11;
default.dec_levels = 7;
default.wname = 'db25';
default.sigma = 2.4;
default.rot_axis_offset = [];
default.rot_axis_tilt = 0.001;
default.parfolder = '';
default.write_to_scratch = 0;
default.write_flatcor = 0;
default.write_phase_map = 0; 
default.write_sino = 0; 
default.write_sino_phase = 0; 
default.write_reco = 1; 
default.write_float = 0; 
default.write_float_binned = 1; 
default.write_8bit = 0;
default.write_8bit_binned = 1;
default.write_16bit = 0; 
default.subfolder_reco = '';
default.gpu_ind = 1;

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2016 commissioning

raw = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/';
default.raw_roi = [141 1940];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg5gd_02_1w'];para(nn).rot_axis_offset = 2.0;
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg10gd_38_1w'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg10gd_41_2w'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg10gd_44_3w'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg10gd_50_4w']; % rec o problems CHECK
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg5gd_02_1w'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg5gd_13_2w'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg5gd_21_3w'];
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw, 'mg5gd_25_4w'];


nn = nn + 1;para(nn) = default;para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160915_000_synload/raw/mg10gd_50_4w';
para(nn).raw_roi = [121 2240];

'/asap3/petra3/gpfs/p05/2016/commissioning/c20160920_000_diana/raw/Mg-10Gd39_1w';
'/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/mg5gd_21_3w';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( nums, doreco, print_field, para)
