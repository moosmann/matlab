function [vol, reco_metric] = find_rot_axis_offset( tomo, proj)
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
% Written by Julian Moosmann. Last modification: 2018-04-11
%
% [vol, reco_metric] = find_rot_axis_offset( tomo, proj );

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slice = assign_from_struct( tomo, 'slice', [] );
vol_shape = assign_from_struct( tomo, 'vol_shape', [] );
vol_size = assign_from_struct( tomo, 'vol_size', [] );
offset = double( assign_from_struct( tomo, 'offset', -10:10 ));
offset_shift = assign_from_struct( tomo, 'offset_shift', 0 );
tilt = assign_from_struct( tomo, 'tilt', 0 );
lamino = assign_from_struct( tomo, 'lamino', 0 );
fixed_tilt = assign_from_struct( tomo, 'fixed_tilt', 0 );
take_neg_log = assign_from_struct( tomo, 'take_neg_log', 1 );
number_of_stds = assign_from_struct( tomo, 'number_of_stds', 4 );
butterworth_filtering = assign_from_struct( tomo.butterworth_filter, 'apply', 0 );
vert_shift = assign_from_struct( tomo, 'vert_shift', 0 );
mask_rad = 0.95;
mask_val = 0;
filter_histo_roi = 0.25;

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
vol = zeros(vol_shape(1), vol_shape(2), numel(offset));
for nn = 1:numel(reco_metric)
    reco_metric(nn).val = zeros( numel(offset), 1);
end

if ~lamino
    tomo.tilt_camera = tilt;
    tomo.tilt_lamino = fixed_tilt;
else
    tomo.tilt_camera = fixed_tilt;
    tomo.tilt_lamino = tilt;
end

% Backprojection
for nn = 1:numel( offset )
    tomo.rot_axis.offset = offset(nn) + offset_shift + eps;
    
    %% Reco
    switch lower( tomo.reco_mode )
        case '3d'
            im = astra_parallel3D( tomo, permute( sino, [1 3 2]) );
        case 'slice'
            im = astra_parallel2D( tomo, permute( sino, [3 1 2]) );
    end
    vol(:,:,nn) = FilterHisto(im, number_of_stds, filter_histo_roi);
    %vol(:,:,nn) = FilterOutlier( im, 0.01 );
    
    %% Metrics    
    im = double( MaskingDisc( im, mask_rad, mask_val) ) * 2^16;
    % mean    
    reco_metric(1).val(nn) = mean2( im );
    % mean abs
    reco_metric(2).val(nn) = mean2( abs( im ) );
    % mean negativity
    reco_metric(3).val(nn) = - mean( im( im <= 0 ) );
    % isotropic gradient
    [g1, g2] = gradient(im);
    reco_metric(4).val(nn) = mean2( sqrt( g1.^2 + g2.^2 ) );
    % laplacian
    reco_metric(5).val(nn) = mean2( abs( del2( im ) ) );
    % entropy
    p = histcounts( im(:) );
    p = p(p>0);
    p = p / sum( p );
    reco_metric(6).val(nn) = -sum( p .* log2( p ) );
    % entropy built-in
    reco_metric(7).val(nn) = -entropy( im );

end

% Normalize for ease of plotting and comparison
for nn = 1:numel(reco_metric)
    reco_metric(nn).val = normat( reco_metric(nn).val );
end
