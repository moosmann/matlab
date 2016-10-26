% Tomographic reconstruction of Shepp-Logan phantom using filtered
% backprojection. Check FBP filter when using ASTRA

clear all
aclear

%% Parameter
N = 128;
M = N +0;
vol_shape_3d = [N, M, 1];
vol_shape_geom = [vol_shape_3d(2), vol_shape_3d(1), vol_shape_3d(3)];

num_angles = 1 * N;
%theta = 0:pi/NumProj:pi-1/NumProj;
angles = linspace2(0, 1*pi, num_angles); % excludes upper limit

det_row_count = 1;
det_col_count = N + 100; % for 3D use square detector
det_width = 1 ;% for 3D use quadratic detector

%% Phantom and projections
P = zeros(N,M);
x = ceil(N/4):floor(3/4*N);
y = ceil(M/4):floor(3/4*M);
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
sino = sino + 0.01 * ( max(sino(:)) - min(sino(:)) ) * ( rand( size(sino) ) - 0.5 ) ;

%% Reco
rotation_axis_offset = 0;
vol_shape = vol_shape_3d;
vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];
pixel_size = [det_width, det_width];
link_data = 0;
vol_shape = vol_shape .* [1 1 1];
% unfiltered
bp = astra_parallel3D(sino, angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

% filtered
% Design filter
filt_len = size(sino,1);
filt = iradonDesignFilter('Ram-Lak', filt_len, 1);

% Fourier transform, filter mutliplication, inverse Fourier transform
sinof = real( ifft( bsxfun(@times, fft( sino, [], 1), filt), [], 1, 'symmetric') );
fbp = astra_parallel3D(sinof, angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

%% Butterworth
[b, a] = butter(1, 0.5);
filt_len_bw = 2 * filt_len;
bw = freqz(b, a, filt_len_bw );
filt_bw = iradonDesignFilter('Ram-Lak', filt_len_bw, 1);
filt_bw = filt_bw .* bw;
pad_vec = zeros([1, ndims(sino)] );
pad_vec(1) = filt_len_bw - filt_len;
sinop = padarray(sino, pad_vec,  0, 'post');
sinof_bw = real( ifft( bsxfun(@times, fft( sinop, [], 1), filt_bw), [], 1, 'symmetric') );
sinof_bw = sinof_bw(1:filt_len,:);
fbp_bw = astra_parallel3D(sinof_bw, angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

%% idl
s = size(sino);
padding = floor( s(1) * sqrt(2) - s(1) );
sinop = zeros( s(1) + 2 * padding - 0, s(2));

[b, a] = butter(1, 0.5);
bw = freqz(b, a, s(1) );
sino_bw = real( ifft( bsxfun(@times, fft( sino, [], 1), bw), [], 1, 'symmetric') );

sinop(1 + padding:end-padding, :) = sino_bw;

% filtered
filt = iradonDesignFilter('Ram-Lak', size(sinop,1), 1);

% Fourier transform, filter mutliplication, inverse Fourier transform
sinof2 = real( ifft( bsxfun(@times, fft( sinop, [], 1), filt), [], 1, 'symmetric') );
fbp2 = astra_parallel3D(sinof2, angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

%% Figures
show_range = [0 1.1];
%figure('Name','Phantom: Original, FBP, BP')
subplot(2,4,1), imshow(P, show_range), title('Original: P')
subplot(3,4,2), imshow(sino, []), title('sino')
%subplot(2,4,3), imshow(bp, []), title('bp')
subplot(2,4,3), imshow(fbp_bw -fbp, []), title('fbp bw -fbp')
subplot(2,4,4), imshow(fbp2, show_range), title('fbp idl')
subplot(2,4,5), imshow(fbp, show_range), title('fbp')
subplot(2,4,6), imshow(fbp_bw, show_range), title('fbp bw')
show_range_err = [0 0.5];
subplot(2,4,7), imshow(abs(fbp-P), show_range_err), title('fbp - P')
subplot(2,4,8), imshow(abs(fbp_bw-P), show_range_err), title('fbp bw - P')

%imshow( fbp - P, [], 'InitialMagnification', 'fit' )
%title(sprintf('Min: %8.4f, Max: %8.4f', min(fbp(:)), max(fbp(:)) ))

domain(fbp)
domain(fbp_bw)
domain(fbp2)
%domain(fbp_bw -fbp, 1, 'fbp bw - fbp')


% [R, Xp] = radon(P,theta);
% RotAxis = ceil(size(R,1)/2);
% 
% %% Inverse radon transformation
% freqScal = 1;
% OutputSize = size(R,1);
% [I1,h] = iradon(R,theta,'linear','Ram-Lak',freqScal,OutputSize);
% 
% %% Plotting
% %figure('Name','Sinogram')
% %imshow(R,[])
% figure('Name','Phantom: Original, FBP, BP')
% subplot(1,3,1), imshow(P), title('Original')
% subplot(1,3,2), imshow(I1), title('FBP')
