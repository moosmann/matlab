function p05_reco_test_ringfilter(nums, doreco, print_field)

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
default.fbp_filter_padding = 1;
default.fbp_filter_type = 'Ram-Lak';'linear';
default.fpb_freq_cutoff = 1; 
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
default.write_reco = 1; 
default.write_float = 1; 
default.write_float_binned = 1; 
default.write_8bit = 0;
default.write_8bit_binned = 0;
default.write_16bit = 0; 
default.subfolder_reco = '';
default.gpu_ind = 1;

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2016 nov
raw = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/';

nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'mah_01'];
para(nn).raw_roi = [293 1492];
para(nn).rot_axis_offset = -135.5 / para(nn).raw_bin;
para(nn).rot_axis_tilt = -0.003;
para(nn).write_sino = 1;
para(nn).ring_filter_method = 'wavelet-fft';
para(nn).parfolder = ['test_ringfilt_' para(nn).ring_filter_method];

nn = nn + 1;para(nn) = para(nn-1);
para(nn).dec_levels = 6;
para(nn).wname = 'db30';
para(nn).parfolder = sprintf( 'test_ringfilt_%s_%s_decNum%u', para(nn).ring_filter_method, para(nn).wname, para(nn).dec_levels );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).ring_filter_method = 'jm';
para(nn).parfolder = ['test_ringfilt_' para(nn).ring_filter_method];

nn = nn + 1;para(nn) = para(nn-1);
para(nn).ring_filter = 0;        
para(nn).ring_filter_method = 'none';
para(nn).parfolder = ['test_ringfilt_' para(nn).ring_filter_method];


%% 2017 may
raw = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';
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
para(nn).ring_filter = 1;
%para(nn).raw_roi = [141 1940];
para(nn).raw_roi = [701 1100];
para(nn).parfolder = 'test_ringfilt_jm';

nn = nn + 1;para(nn) = para(nn-1);
para(nn).rot_axis_offset = 2.0;
para(nn).wname = 'db15';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:3;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).wname = 'db15';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:7;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).rot_axis_offset = 2.0;
para(nn).wname = 'db25';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:3;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).wname = 'db25';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:7;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).rot_axis_offset = 2.0;
para(nn).wname = 'db30';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:3;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).wname = 'db30';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:7;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).rot_axis_offset = 2.0;
para(nn).wname = 'bior2.8';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:3;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).wname = 'bior2.8';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:7;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', para(nn).wname, sprintf('%u', para(nn).dec_levels) );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).wname = 'db30';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:7;
para(nn).fbp_filter_type = 'linear';
para(nn).fpb_freq_cutoff = 0.5;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s_fbpFilt%s_fbpcut%02u', para(nn).wname, sprintf('%u', para(nn).dec_levels), para(nn).fbp_filter_type, para(nn).fpb_freq_cutoff*10 );

nn = nn + 1;para(nn) = para(nn-1);
para(nn).wname = 'db30';
para(nn).sigma = 2.4;
para(nn).dec_levels = 1:7;
para(nn).fbp_filter_type = 'linear';
para(nn).fpb_freq_cutoff = 1;
para(nn).parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s_fbpFilt%s_fbpcut%02u', para(nn).wname, sprintf('%u', para(nn).dec_levels), para(nn).fbp_filter_type, para(nn).fpb_freq_cutoff*10 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( nums, doreco, print_field, para)
