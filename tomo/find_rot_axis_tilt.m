function [vol, m1, m2, m3, m4, m5] = find_rot_axis_tilt(proj, angles, slice, offset, tilts, take_neg_log, number_of_stds, vol_shape)
% Reconstruct slices from a single sinogram using a range of rotation axis
% tilts.
%
% RETURN
% vol : 3D array. stack of slices with different rotation axis tilts
% m1 : scalar. mean values
% m2 : scalar. mean absolute value
% m3 : scalar. mean non-negative values
% m4 : scalar. mean of isotropic modulus of gradient
% m5 : scalar. mean of Laplacian
% 
% Written by Julian Moosmann. Last modification: 2017-04-20
%
% [vol, m1, m2, m3, m4, m5] = find_rot_axis_tilt(proj, angles, slice, offset, tilts, take_neg_log, number_of_stds, vol_shape)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    slice = [];
end
if nargin < 4
    offset = 0;
end
if nargin < 5
    tilts = -0.005:0.001:0.005;
end
if nargin < 6
   take_neg_log = 1;
end
if nargin < 7
    number_of_stds = 4;
end
if nargin < 8
    vol_shape = [];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_pix, num_row, num_slices] = size( proj );
if isempty( vol_shape )
    vol_shape = [num_pix, num_pix, 1];
else
    vol_shape(3) = 1;
end
vol_size = [-num_pix/2 num_pix/2 -num_pix/2 num_pix/2 -0.5 0.5];
astra_pixel_size = 1;
link_data = 1;
roi = 0.25;
if isempty( slice )
    slice = round( num_row / 2 );
end

% Calculate required slab
rot_axis_pos = offset + num_pix / 2;
l = max( rot_axis_pos, abs( num_pix - rot_axis_pos ));
dz = ceil( sin( max( abs( tilts ) ) ) * l ); % maximum distance between sino plane and reco plane
if slice - dz < 0 || slice + dz > num_slices
    fprintf( '\nWARNING: Inclination of reconstruciont plane exceeds sinogram volume. Better choose a more central slice or a smaller tilts.')
end

% Slab
y_range = slice + (-dz:dz);
sino = proj(:, y_range, :);

% Ramp filter
filt = iradonDesignFilter('Ram-Lak', 2 * num_pix, 1);

% Butterworth filter
[b, a] = butter(1, 0.5);
bw = freqz(b, a, numel(filt) );
filt = filt .* bw;

% Apply filters
% Apply filters
sino = padarray( NegLog(sino, take_neg_log), [num_pix 0 0 ], 'symmetric', 'post');
sino = real( ifft( bsxfun(@times, fft( sino, [], 1), filt), [], 1, 'symmetric') );
sino = sino(1:num_pix,:,:);

% Preallocation
vol = zeros(num_pix, num_pix, numel(tilts));
m1 = zeros(1, numel(tilts));
m2 = m1;
m3 = m1;
m4 = m1;
m5 = m1;
roi1 = IndexParameterToRange( 0, num_pix );

% fprintf('y range: ')
% fprintf(' %u', y_range)
% fprintf(' dz : %g', dz )
% fprintf( 'angle : %g rad, %g degree', max(abs(tilts)), max( abs(tilts) ) * 180 / pi)

% Backprojection
for nn = 1:numel( tilts )
    tilt = tilts(nn);
    
    % Reco
    im = FilterHisto(astra_parallel3D( permute( sino, [1 3 2]), angles, offset, vol_shape, vol_size, astra_pixel_size, link_data, tilt), number_of_stds, roi);    
    vol(:,:,nn) = im;
    
    % Metrics on ROI
    im_roi = im( roi1, roi1);
    % Mean
    m1(nn) = mean2( im_roi );
    % Mean abs
    m2(nn) = mean2( abs( im_roi ) );
    % Mean negativity
    m3(nn) = - mean( im_roi( im_roi <= 0 ) );
    % isotropic gradient
    [g1, g2] = gradient(im_roi);
    m4(nn) = mean2( sqrt( g1.^2 + g2.^2 ) );
    % laplacian
    m5(nn) = mean2( del2( im_roi ) );

end

%m1 = m1 ./ mean2( sino );
m1 = normat( m1 );
m2 = normat( m2 );
m3 = normat( m3 );
m4 = normat( m4 );
m5 = normat( m5 );
