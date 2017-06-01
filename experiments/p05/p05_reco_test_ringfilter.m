function p05_reco_test_ringfilter(nums, doreco, print_field)

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
nn = 0;
raw = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';
default.scan_path = '';
default.raw_roi = [299 1492];
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
default.fbp_filter_padding = 0;
default.ring_filter = 1;
default.ring_filter_method = 'wavelet-fft';
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
default.write_reco = 0; 
default.write_float = 1; 
default.write_float_binned = 1; 
default.write_8bit = 0;
default.write_8bit_binned = 0;
default.write_16bit = 0; 
default.subfolder_reco = '';
default.gpu_ind = 1;

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Combined wavelet FFT
nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_00'];
para(nn).ring_filter_method = 'wavelet-fft';
para(nn).wname = 'db25';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1); 
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', para(nn).dec_levels);

% Simple JM
nn = nn + 1;para(nn) = para(nn-1);
para(nn).ring_filter_method = 'jm';
para(nn).parfolder = 'test_ringfilt_jm';

% None
nn = nn + 1;para(nn) = para(nn-1);
para(nn).ring_filter = 0;
para(nn).parfolder = 'test_ringfilt_none';

% DB30
nn = nn + 1;para(nn) = default; 
para(nn).ring_filter_method = 'wavelet-fft';
para(nn).wname = 'db30';
para(nn).sigma = 2.4;
para(nn).dec_levels = 3;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', para(nn).wname, para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', para(nn).wname, para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', para(nn).wname, para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', para(nn).wname, para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = para(nn-1).dec_levels + 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', para(nn).wname, para(nn).dec_levels);

%% Synload
nn = nn + 1;para(nn) = para(nn-1);para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/mg5gd_02_1w';
para(nn).rot_axis_offset = 2.0;
para(nn).ring_filter = 1;
para(nn).wname = 'db30';
para(nn).sigma = 2.4;
para(nn).dec_levels = 5;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', para(nn).wname, para(nn).dec_levels);

nn = nn + 1;para(nn) = para(nn-1);para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/mg5gd_02_1w';
para(nn).ring_filter = 0;
para(nn).parfolder = 'test_ringfilt_jm';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( nums, doreco, print_field, para)
