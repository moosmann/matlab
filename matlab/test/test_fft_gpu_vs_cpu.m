fprintf( '\nFTT TEST' )
clear all

warning( 'off', 'parallel:gpu:device:DeviceDeprecated' )
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, use to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D';'slice';  % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = [-0.5   1 -0.5 0.5 -0.5 0.5];%[];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_full_range = [];% 2 * pi ;[]; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle_offset = pi; % global rotation of reconstructed volume
tomo.rot_axis_offset = []; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.flip_scan_position = 0; % for debugging
tomo.rot_axis_tilt_camera = 0; % in rad. camera tilt w.r.t rotation axis.
tomo.rot_axis_tilt_lamino = 0; % in rad. lamino tilt w.r.t beam.
tomo.rot_axis_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
tomo.rot_axis_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
tomo.fbp_filter_type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
tomo.fbp_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
tomo.fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
tomo.fbp_filter_padding_method = 'symmetric';
tomo.butterworth_filter = 0; % use butterworth filter in addition to FBP filter
tomo.butterworth_filter_order = 1;
tomo.butterworth_filter_frequ_cutoff = 0.9;
tomo.astra_pixel_size = 1; % detector pixel size for reconstruction: if different from one 'tomo.vol_size' must to be ajusted, too!
tomo.take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
tomo.algorithm = 'fbp';'sirt'; 'cgls';'sart';'em';'fbp-astra'; % SART/EM only work for 3D reco mode
tomo.iterations = 40; % for iterateive algorithms: 'sirt', 'cgls', 'sart', 'em'
tomo.MinConstraint = []; % sirt3D/sirt2d/sart2d only. If specified, all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions, for example.
tomo.MaxConstraint = [];

%raw_path = '/asap3/petra3/gpfs/p07/2020/data/11010172/raw/swerim_21_12_oh_a/tiff0000/';
fn = '/asap3/petra3/gpfs/p07/2020/data/11010172/raw/swerim_21_12_oh_a/tiff0000/swerim_21_12_oh_a_0001045_img.tif';
warning( 'off', 'imageio:tifftagsread:expectedTagDataFormat'  )

if ~exist( 'im0', 'var' )
    imu = imread( fn, 'tif' );
    ims = single( imu' );
end
im_size = size( ims );
proj_shape1 = im_size(1);
cl0 = class( imu );
cl = class( ims );
fprintf( '\n size : %u %u', im_size )
fprintf( '\n class : %s, %s', cl0, cl )


filt = iradonDesignFilter( tomo.fbp_filter_type, (1 + tomo.fbp_filter_padding) * proj_shape1, tomo.fbp_filter_freq_cutoff);
take_neg_log = tomo.take_neg_log;
padding = tomo.fbp_filter_padding;
padding_method = tomo.fbp_filter_padding_method;

if 0
    g = gpuDevice;
    g.reset
    im = ims;
    vol = cat( 3, im, im, im, im, im );
end

% CPU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    fprintf( '\nCPU ' )
    im = ims;
    
    tic
    im = padarray( im, padding * [proj_shape1 0 0], padding_method, 'post' );
    im = fft( im, [], 1);
    im = bsxfun( @times, im, filt );
    im = real( ifft( im, [], 1, 'symmetric') );
    im = im(1:proj_shape1,:,:);
    tcpu = toc;
    fprintf( '%f s', tcpu )
end

% GPU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    fprintf( '\nGPU ' )
    
    tic;
    im = gpuArray( ims );
    im = padarray( im, padding * [proj_shape1 0 0], padding_method, 'post' );
    im = fft( im, [], 1);
    im = bsxfun( @times, im, filt );
    im = real( ifft( im, [], 1, 'symmetric') );
    im = im(1:proj_shape1,:,:);
    img = gather( im );
    tgpu = toc;
    
    fprintf( '%f s', tgpu )
end

if 1
    test_gpu_slab_single_var( ims )
end

if 1
    test_gpu_slab_two_var( ims )
end

fprintf( '\n' )

function avail = gpuavail( what, start)
g = gpuDevice;
tot = g.TotalMemory;
avail = g.AvailableMemory;
if nargin < 2
    start = tot;
end
fprintf( '\n %8s %8.0f MiB, %8.0f MiB, %8.0f MiB', what, avail / 1024^2, ( tot - avail ) / 1024^2, ( start - avail ) / 1024^2 )
end

function test_gpu_slab_single_var(im)

fprintf( '\n\nGPU slab single var ' )

padding = 1;padding_method = 'symmetric';

vol = cat( 3, im, im, im, im, im );
vol_mem =  numel( vol ) * 4 / 1024^2;

fprintf( '\n vol mem : %.0f MiB', vol_mem )
g = gpuDevice(2);
fprintf( '\n total %.0f MiB', g.TotalMemory / 1024^2 )
fprintf( '\n avail %.0f MiB', g.AvailableMemory / 1024^2 )

filt = iradonDesignFilter( 'Ram-Lak', 2*size(im,1),1);
filtg = gpuArray( filt );

tic;
s = gpuavail( 'init' );
gpuavail( 'start', s );

volg = gpuArray( vol );
gpuavail( 'gpuArray', s );

volg = padarray( volg, padding * [size( im, 1) 0 0], padding_method, 'post' );
gpuavail( 'pad', s );

volg = fft( volg, [], 1);
gpuavail( 'fft', s );

volg = bsxfun( @times, volg, filtg );
gpuavail( 'filt', s);

volg = ifft( volg, [], 1, 'symmetric');
gpuavail( 'ifft', s );

volg = real( volg );
gpuavail( 'real', s );

volg = volg(1:size( im, 1),:,:);
gpuavail( 'crop', s );

vol2 = gather( volg );
gpuavail( 'gather', s );

tgpu = toc;
fprintf( '\n' )
fprintf( '%f s, %f s', tgpu, tgpu / size( vol, 3 ) )
end

function test_gpu_slab_two_var( im )

fprintf( '\n\nGPU slab two var' )

vol = cat( 3, im, im, im, im, im );

padding = 1;padding_method = 'symmetric';

g = gpuDevice(2);
vol_mem =  numel( vol ) * 4 / 1024^2;
fprintf( '\n vol mem : %.0f MiB', vol_mem )
fprintf( '\n total %.0f MiB', g.TotalMemory / 1024^2 )
fprintf( '\n avail %.0f MiB', g.AvailableMemory / 1024^2 )
filt = iradonDesignFilter( 'Ram-Lak', 2*size(im,1),1);
filtg = gpuArray( filt );
tic;
s = gpuavail( 'init' );
gpuavail( 'start', s );

volg = gpuArray( vol );
gpuavail( 'gpuArray', s );

volg2 = padarray( volg, padding * [size( im, 1) 0 0], padding_method, 'post' );
gpuavail( 'pad', s );

volg = fft( volg2, [], 1);
gpuavail( 'fft', s );

volg2 = bsxfun( @times, volg, filtg );
gpuavail( 'filt', s);

volg = ifft( volg2, [], 1, 'symmetric');
gpuavail( 'ifft', s );

volg2 = real( volg );
gpuavail( 'real', s );

volg = volg2(1:size( im, 1),:,:);
gpuavail( 'crop', s );

vol2 = gather( volg );


tgpu = toc;

fprintf( '%f s, %f s', tgpu, tgpu / size( vol, 3 ) )
end