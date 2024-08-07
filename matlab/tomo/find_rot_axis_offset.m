function [vol,reco_metric] = find_rot_axis_offset(tomo,proj,par)
% Reconstruct slices from sinogram for a range of rotation axis position
% offsets.
%
% RETURN
% vol : 3D array. stack of slices with different rotation axis position
% offsets
% reco_metric : struct containing different metrics: mean of all values,
% mean of all absolute values, mean non-negative values, mean of isotropic
% modulus of gradient, mean of Laplacian, entropy
%
% Written by Julian Moosmann.
%
%[vol,reco_metric] = find_rot_axis_offset(tomo,proj,par)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slice = assign_from_struct( tomo, 'slice', [] );
vol_shape = assign_from_struct( tomo, 'vol_shape', [] );
vol_size = assign_from_struct( tomo, 'vol_size', [] );
offset = double( assign_from_struct( tomo, 'offset', -10:10 ));
%offset_shift = assign_from_struct( tomo, 'offset_shift', 0 );
tilt = assign_from_struct( tomo, 'tilt', 0 );
lamino = assign_from_struct( tomo, 'lamino', 0 );
fixed_tilt = assign_from_struct( tomo, 'fixed_tilt', 0 );
take_neg_log = assign_from_struct( tomo, 'take_neg_log', 1 );
butterworth_filtering = assign_from_struct( tomo, 'butterworth_filter', 0 );
vert_shift = assign_from_struct( tomo, 'vert_shift', 0 );
mask_rad = 0.95;
mask_val = 0;
%number_of_stds = assign_from_struct( tomo, 'number_of_stds', 4 );
%filter_histo_roi = 0.25;
rect = assign_from_struct(tomo, 'rot_axis_offset_metric_roi', []);               

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_pix, num_row, ~] = size( proj );
if isempty( vol_shape )
    vol_shape = [num_pix, num_pix, 1];
else
    vol_shape(3) = 1;
end
tomo.vol_shape = vol_shape;
if isempty( vol_size )
    vol_size = [-num_pix/2 num_pix/2 -num_pix/2 num_pix/2 -0.5 0.5];
else
    vol_size(5) = -0.5;
    vol_size(6) = 0.5;
end
tomo.vol_size = vol_size;
if isempty( slice )
    slice = round( num_row / 2 );
else
    if slice < 1 && slice >= 0
        slice = floor((size( proj, 2 ) - 1) * slice + 1 );
    end
end
tomo.slice = slice;

% Calculate required slab size: tilt condition
rot_axis_pos = offset + vol_shape(1) / 2;
l = max( max( rot_axis_pos) , max( abs( vol_shape(1) - rot_axis_pos ) ) );
% maximum distance between sino plane and reco plane
dz = ceil( sin( max( abs( [ tilt, fixed_tilt] ) ) ) * l );
if slice - dz < 0 || slice + dz > num_row
    fprintf( '\nWARNING: Inclination of reconstruction plane, slice %u, exceeds sinogram volume. Better choose a more central slice or a smaller tilts.', slice)
end
% Calculate required slab size: Spiral CT condition
if numel( vert_shift ) > 1
    dz = dz + floor( max( abs( tomo.vert_shift ) ) );
end
if slice - dz < 0 || slice + dz > num_row
    fprintf( '\nWARNING: Spiral CT requires larger sinogram volume. Better choose a more central slice or a smaller tilts.')
end

% Slab
slice_range = slice + (-dz:dz);
sino = proj(:, slice_range, :);
%fprintf( '\n sino slab size : %u %u %u', size(sino) )
fprintf(' \n Filter sino')
if strcmpi( tomo.algorithm, 'fbp' )
    % Ramp filter
    filt = iradonDesignFilter('Ram-Lak', 2 * num_pix, 0.9);
    
    % Butterworth filter
    if butterworth_filtering
        [b, a] = butter(1, 0.5);
        bw = freqz(b, a, numel(filt) );
        filt = filt .* bw;
    end
    
    % Apply filters
    sino = padarray( NegLog(sino, take_neg_log), [num_pix 0 0 ], 'symmetric', 'post');
    sino = real( ifft( bsxfun(@times, fft( sino, [], 1), filt), [], 1, 'symmetric') );
    sino = sino(1:num_pix,:,:);
end

% Metrics
reco_metric(1).name = 'mean';
reco_metric(2).name = 'abs';
reco_metric(3).name = 'neg';
reco_metric(4).name = 'iso-grad';
reco_metric(5).name = 'laplacian';
reco_metric(6).name = 'entropy';
reco_metric(7).name = 'entropy-ML';

% Preallocation
vol = zeros(vol_shape(1), vol_shape(2), numel(offset), 'single');
for nn = 1:numel(reco_metric)
    reco_metric(nn).val = zeros( numel(offset), 1);
end

if ~lamino
    tomo.rot_axis_tilt_camera = tilt;
    tomo.rot_axis_tilt_lamino = fixed_tilt;
else
    tomo.rot_axis_tilt_camera = fixed_tilt;
    tomo.rot_axis_tilt_lamino = tilt;
end

% Permute sino for ASTRA
tomo_reco_mode = lower( tomo.reco_mode );
switch tomo_reco_mode
    case '3d'
        sino = permute( sino, [1 3 2]);
    case 'slice'
        sino = permute( sino, [3 1 2]);
end

%% Reco loop over different offset
reco_metric_mean = zeros( numel(offset), 1);
reco_metric_abs = zeros( numel(offset), 1);
reco_metric_neg = zeros( numel(offset), 1);
reco_metric_grad = zeros( numel(offset), 1);
reco_metric_lap = zeros( numel(offset), 1);
reco_metric_ent = zeros( numel(offset), 1);
reco_metric_entml = zeros( numel(offset), 1);
%rect(2) + (1:rect(4)), rect(1) + (1:rect(3))

% num_gpu = numel(tomo.gpu_index);
% if isempty(num_gpu)
%     num_gpu = gpuDeviceCount;
% end
% gpu_index = tomo.astra_gpu_index;
% fprintf('\n ASTRA GPU index (count from 0): ')
% fprintf(' %u,', gpu_index)

gpu_index_list = p_gpu_list(sino, tomo, par, 1);
num_gpu_used = numel(gpu_index_list);

% nn-dim struct required for parfor loop
%tomo_par = struct([]);
tomo.astra_gpu_index = [];
for nn = numel(offset):-1:1
    tomo_par(nn) = tomo;
    tomo_par(nn).rot_axis_offset = offset(nn);
    mm = mod(nn - 1, num_gpu_used) + 1;
    tomo_par(nn).astra_gpu_index = gpu_index_list(mm);
    %fprintf('\n  par pool index: %2u, gpu list index: %2u, gpu index: %u', nn, mm, gpu_index_list(mm))
end

fprintf('\n Offset par loop:\n')
fprintf('\n slice: %u',slice)
if isempty(rect)
    fprintf('\n roi: []')
else
    fprintf('\n roi: ')
    fprintf(' %u',rect)
end
parfor (nn = 1:numel(offset), 2*num_gpu_used)
%parfor (nn = 1:numel(offset), 2 * num_gpu)
%for nn = 1:numel(offset)
    
    % Reco
    tomo_par_nn = tomo_par(nn);
    %     mm = mod(nn, num_gpu) + 1;
%     tomo_par_nn.astra_gpu_index = gpu_index_list(mm); %#ok<PFBNS>
%     fprintf('\n  par index: %2u, gpu list index: %2u, gpu index: %u', nn, mm, gpu_index_list(mm))
    im = [];
    switch tomo_reco_mode
        case '3d'
            im = astra_parallel3D( tomo_par_nn, sino );
        case 'slice'
            gpu_index = mod( nn,  num_gpu_used ) + 1;            
            im = astra_parallel2D( tomo_par_nn, sino, gpu_index);
            %im = astra_parallel2D( tomo_par_nn, sino);
    end
    vol(:,:,nn) = im;
    %vol(:,:,nn) = FilterHisto(im, number_of_stds, filter_histo_roi);
    %vol(:,:,nn) = FilterOutlier( im, 0.01 );
    
    % ROI
    if isempty(rect)
        im = MaskingDisc( im, mask_rad, mask_val);
    else
        im = im(rect(2) + (1:rect(4)), rect(1) + (1:rect(3)));
    end
    im = 2^16 * double( im );
    
    % Metrics
    % mean
    reco_metric_mean(nn) = mean2( im );
    % mean abs
    reco_metric_abs(nn) = mean2( abs( im ) );
    % mean negativity
    reco_metric_neg(nn) = - mean( im( im <= 0 ) );
    % isotropic gradient
    [g1, g2] = gradient(im);
    reco_metric_grad(nn) = mean2( sqrt( g1.^2 + g2.^2 ) );
    % laplacian
    reco_metric_lap(nn) = mean2( abs( del2( im ) ) );
    % entropy
    p = histcounts( im(:) );
    p = p(p>0);
    p = p / sum( p );
    reco_metric_ent(nn) = -sum( p .* log2( p ) );
    % entropy built-in
    reco_metric_entml(nn) = -entropy( im );
end
aclear
reco_metric(1).val = reco_metric_mean;
reco_metric(2).val = reco_metric_abs;
reco_metric(3).val = reco_metric_neg;
reco_metric(4).val = reco_metric_grad;
reco_metric(5).val = reco_metric_lap;
reco_metric(6).val = reco_metric_ent;
reco_metric(7).val = reco_metric_entml;

% Normalize for ease of plotting and comparison
for nn = 1:numel(reco_metric)
    reco_metric(nn).val = normat( reco_metric(nn).val );
end
