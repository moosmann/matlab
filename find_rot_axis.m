function [vol, sino, m1, m2, m3, m4, m5, com] = find_rot_axis(proj, angles, offsets, slice, take_neg_log)

%% Default arguments
if nargin < 4
    slice = [];
end
if nargin < 5
    take_neg_log = 0;
end
num_pix = size(proj, 1);
subvol_shape = [num_pix, num_pix, 1];
subvol_size = [-num_pix/2 num_pix/2 -num_pix/2 num_pix/2 -0.5 0.5];
astra_pixel_size = 1;
link_data = 1;
rot_axis_tilt = 0;
number_of_stds = 3;
roi = 0.25;

%% Main
if isempty( slice )
    slice = round( size( proj, 2) / 2 );
end

sino = squeeze( proj(:, slice, :));
[~, x] = meshgrid( 1:size(sino,2), 1:size(sino,1));
com = mean2( x.*sino );

% Ramp filter
filt = iradonDesignFilter('Ram-Lak', num_pix, 0.9);

% Butterworth filter
[b, a] = butter(1, 0.5);
bw = freqz(b, a, numel(filt) );
filt = filt .* bw;

% Apply filters
sino = real( ifft( bsxfun(@times, fft( NegLog(sino, take_neg_log), [], 1), filt), [], 1, 'symmetric') );

% Preallocation
vol = zeros(num_pix, num_pix, numel(offsets));
m1 = zeros(1, numel(offsets));
m2 = m1;
m3 = m1;
m4 = m1;
m5 = m1;
roi1 = IndexParameterToRange( 0, num_pix );

% Backprojection
for nn = 1:numel( offsets )
    offset = offsets(nn);
    
    % Reco
    im = FilterHisto(astra_parallel3D( sino, angles, offset, subvol_shape, subvol_size, astra_pixel_size, link_data, rot_axis_tilt), number_of_stds, roi);    
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