% Tomographic reconstruction of Shepp-Logan phantom using filtered
% backprojection. Check FBP filter when using ASTRA

clear all
aclear

%% Parameter
N = 1000;
M = N +0;
vol_shape_3d = [N, M, 1];
vol_shape_geom = [vol_shape_3d(2), vol_shape_3d(1), vol_shape_3d(3)];

num_angles = round(1.5 * N);
%theta = 0:pi/NumProj:pi-1/NumProj;
angles = linspace2(0, 1*pi, num_angles); % excludes upper limit

det_row_count = 1;
det_col_count = N + 100; % for 3D use square detector
det_width = 1 ;% for 3D use quadratic detector

%% Phantom and projections
P = zeros(N,M);
x = ceil((N-1)/4+1):floor(3/4*(N-1)+1);
y = ceil((M-1)/4+1):floor(3/4*(M-1)+1);
P(x,y) = 1;
P = P + 0.0;

%% 3D parallel FP

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape_geom );

% projection geometry
det_spacing_x = det_width;
det_spacing_y = det_width;
proj_geom = astra_create_proj_geom('parallel3d', det_spacing_x, ...
    det_spacing_y, det_row_count, det_col_count, angles);

% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, P);

% create sino
sino_id = astra_mex_data3d('create','-sino', proj_geom, 0);
    
% create struct
cfg = astra_struct('FP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;
cfg.option.GPUindex = 1; 

% forward projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
sino = astra_mex_data3d('get',sino_id);

%% Noise
sino0 = sino;
%sino = sino + 0.01 * ( max(sino(:)) - min(sino(:)) ) * ( rand( size(sino) ) - 0.5 ) ;

%% Reco %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotation_axis_offset = 0;
vol_shape = vol_shape_3d;
vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];
pixel_size = [det_width, det_width];
link_data = 0;
vol_shape = vol_shape .* [1 1 1];
freq_cutoff = 1;
pad = @(sino) padarray(padarray(sino, [2*ceil((N*sqrt(2) - size(sino,1))) 0],0, 'pre'), [2*ceil((N*sqrt(2) - size(sino,1))) 0],0, 'post');
pad = @(sino) padarray(padarray(sino, [ceil((2*N - size(sino,1))/2) 0],0, 'pre'), [ceil((2*N - size(sino,1))/2) 0],0, 'post');
%pad = @(sino) sino;
%% unfiltered
bp = astra_parallel3D(pad(sino), angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);
len = size(sino,1);

%% Ram-Lak
[b,a]=butter(1,0.5);
%sinop=pad(filter(b,a,pad(sino)));
sinop = sino;
%sinop = padarray( sino, [ 2^nextpow2(len); - len 0],0, 'post');
filt_len = size(sinop,1);
filt_rl = iradonDesignFilter('Ram-Lak', filt_len, freq_cutoff);
sinof_rl = real( ifft( bsxfun(@times, fft( sinop, [], 1), filt_rl), [], 1, 'symmetric') );
sino_rl = sinof_rl(1:len,:);
fbp_rl = astra_parallel3D(pad(sino_rl), angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

%% Ram-Lak plus Butterworth
[b, a] = butter(1, 0.5);
bw = freqz(b, a, filt_len );
filt_bw = filt_rl  .* bw;
sinof_bw = real( ifft( bsxfun(@times, fft( sinop, [], 1), filt_bw), [], 1, 'symmetric') );
sinof_bw = sinof_bw(1:len,:);
fbp_rl_bw = astra_parallel3D(pad(sinof_bw), angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

%% Linear
filt_lin = iradonDesignFilter('linear', filt_len, freq_cutoff);
sinof_lin = real( ifft( bsxfun(@times, fft( sinop, [], 1), filt_lin), [], 1, 'symmetric') );
sinof_lin = sinof_lin(1:len,:);
fbp_lin = astra_parallel3D(pad(sinof_lin), angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

%% idl
s = size(sino);
padding = floor( s(1) * sqrt(2) - s(1) );
sinop2 = zeros( s(1) + 2 * padding - 0, s(2));
[b, a] = butter(1, 0.5);
bw = (freqz(b, a, s(1) ));
sino_bw = real( ifft( bsxfun(@times, fft( sino, [], 1), bw), [], 1, 'symmetric') );
sinop2(1 + padding:end-padding, :) = sino_bw;
%sinop2 = sino_bw;
filt = iradonDesignFilter('Ram-Lak', size(sinop2,1), freq_cutoff);
sino_idl = real( ifft( bsxfun(@times, fft( sinop2, [], 1), filt), [], 1, 'symmetric') );
fbp_idl = astra_parallel3D(sino_idl, angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

%% idl 2
sinop = pad(sino);
s1 = size(sinop,1);
[b, a] = butter(1, 0.5);
bw = freqz(b, a, s1);
filt = iradonDesignFilter('Ram-Lak', s1, freq_cutoff);
sino_idl2 = real( ifft( bsxfun(@times, fft( sinop, [], 1), filt.*bw), [], 1, 'symmetric') );
fbp_idl2 = astra_parallel3D(sino_idl2, angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

%% Figures
show_range = [];[0 1];

figure(1)
subplot(1,3,1), imshow(P, []), title('phantom: P')
subplot(1,3,2), imshow(sino, []), title('sino')
subplot(1,3,3), imshow(bp, []), title('bp')

figure(2)
show_range = [min(min(min([fbp_rl(:), fbp_rl_bw(:), fbp_idl(:), fbp_lin(:)]))), ...
    max(max([fbp_rl(:), fbp_rl_bw(:), fbp_idl(:), fbp_lin(:)]))];
subplot(2,4,1), imshow(fbp_rl, show_range), title('fbp ram-lak')
subplot(2,4,2), imshow(fbp_rl_bw, show_range), title('fbp ram-lak bw')
subplot(2,4,3), imshow(fbp_idl, show_range), title('fbp IDL')
subplot(2,4,4), imshow(fbp_idl2, show_range), title('fbp IDL 2')
%subplot(2,4,4), imshow(fbp_lin, show_range), title('fbp linear')

% Errors
err_rl = abs( fbp_rl - P );
err_rl_bw = abs( fbp_rl_bw - P );
err_idl = abs( fbp_idl - P );
err_idl2 = abs( fbp_idl2 - P );
err_lin = abs( fbp_lin - P );
max_err = max( [ err_rl(:) err_rl_bw(:) err_idl(:) err_lin(:) err_idl2(:)] );
show_range = [0 max(max_err)/2 ];
subplot(2,4,5), imshow(err_rl, show_range), title('fbp ram-lak')
subplot(2,4,6), imshow(err_rl_bw, show_range), title('fbp ram-lak bw')
subplot(2,4,7), imshow(err_idl, show_range), title('fbp IDL')
subplot(2,4,8), imshow(err_idl2, show_range), title('fbp IDL 2')
%subplot(2,4,8), imshow(err_lin, show_range), title('fbp linear')

%% Print
fprintf('\n')
domain(P)
domain(fbp_rl)
domain(fbp_rl_bw)
domain(fbp_idl)
domain(fbp_lin)
% Error
fprintf('                   rl,       rl bw      idl     lin     idl2\n')
fprintf('max errors:   ')
disp( max_err )
fprintf('mean error L1:')
disp( mean( [ err_rl(:) err_rl_bw(:) err_idl(:) err_lin(:) err_idl2(:)] ) );
fprintf('mean error L2:')
disp( sqrt(mean( [ err_rl(:).^2 err_rl_bw(:).^2 err_idl(:).^2 err_lin(:).^2 err_idl2(:).^2  ] ) ) );
domain(err_rl)
domain(err_rl_bw)
domain(err_idl)
domain(err_lin)
