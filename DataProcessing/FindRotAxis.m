function [offset, metric] = FindRotAxis( sino, offset_shift, sampleWidth_to_detectorWidth )
% Determine the center of rotation via Fourier analysis of projection data.
% See articles:
%
% Written by J.Moosmann,2018-08-05. Last version: 2018-08-05

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    %phan = phantom( 128 );
    %sino = radon( phan );
    sino = imread( '/asap3/petra3/gpfs/p05/2018/data/11005553/processed/syn033_68R_Mg10Gd_12w/sino/sino_001000.tif' )';
end
if nargin < 2
    %offset_shift = 0
    load( '/asap3/petra3/gpfs/p05/2018/data/11005553/processed/syn033_68R_Mg10Gd_12w/reco/offset_shift.mat', 'offset_shift' );
end
if nargin < 3
    sampleWidth_to_detectorWidth = 1;
end
if nargin < 4
    offset_init = 0;
end
if nargin < 5
    visual_output = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove lateral sample shift
sino = CropShift( sino, offset_shift );

% Remove stripes
scalfac = mean2( sino ) ./ mean( sino, 1 ); 
sino = bsxfun( @times, sino, scalfac );

% Remove mean
sino = SubtractMean( sino );

[num_pix, num_proj] = size( sino );

% Offset search range
steps = 10;
half_range =  floor( 0.3 * num_pix / 2 );
dx = floor( 2 * half_range / steps );
offset_range = 0:dx:half_range;
offset_range = [-offset_range(end:-1:2) offset_range];
offset_range = offset_init + offset_range;
fprintf( '\n steps : %u', steps )
fprintf( '\n pixels : %u', num_pix )
fprintf( '\n half range : %g', half_range )
fprintf( '\n step size : %g', dx )
fprintf( '\n offset range : ' )
fprintf( ' %g', offset_range )
fprintf( '\n' )

%% Triangular mask
% Fourier coordinates 1D
precision = 'double';
xi  = FrequencyVector( 2 * num_proj, precision, 1 );
eta = FrequencyVector( 2 * num_pix, precision, 1 );
% Fourier coordinates 2D
[xi, eta]  = meshgrid(xi,eta);
mask = eta - abs( sampleWidth_to_detectorWidth * xi ) < 0;
mask = mask & flipud( mask );
sum_mask = sum( mask(:) );

metric = zeros( size( offset_range ) );
if visual_output
    figure( 'Name', 'mask' );
    imsc( fftshift( mask ) )
    axis equal tight
    
    figure( 'Name', 'sino and mask .* FT sino' );
    tmp = log( 10^(-0) + abs( fftshift( fft2( sino ) ) ) );
    mask_plot = 1 + 0.3*(normat( fftshift( mask(1:2:end,1:2:end) ) ) - 0.5);
    imsc( tmp .* mask_plot)
    title( 'sino .* mask. mask scaled for visibility' )
    axis equal tight
    
    drawnow
end

%% Loop over offset range
hsino = figure( 'Name', 'sinogram' );
mask_plot = 1 + 0.3 * (normat( fftshift( mask ) ) - 0.5);
for nn = 1:numel( offset_range )
    offset = offset_range(nn);
    fprintf( '\n step %u, offset: %g', nn, offset )
    if offset > 0
        tmp = cat( 2, sino, cat( 1, sino(end-offset+1:end,:), sino(1:end-offset,:)) );
    else
        tmp = cat( 2, sino, cat( 1, sino(1-offset:end,:), sino(1:-offset,:)) );
    end
    if visual_output
        figure( hsino )
        subplot(2,1,1)
        imsc( tmp )
        title('sinogram')%drawnow;pause(1);
        axis tight equal
    end
    tmp = padarray( tmp, [num_pix, 0], 'symmetric', 'post' );
    tmp = fft2( tmp );
    tmp = abs( tmp );
    
    if visual_output
        subplot(2,1,2)
        tmp_plot = fftshift( log( 1 + abs( tmp ) ) ) .* mask_plot;
        imsc( tmp_plot )
        title( 'mask .* FT sino' )
        axis tight equal
        drawnow;
        pause(1);
    end

    metric(nn) = norm( tmp .* mask, 1 ) / sum_mask;
    fprintf( ', metric: %g ', metric(nn) )
end

%% Fit
offset_range_mean = mean( offset_range );
offset_range_min = min( offset_range );
offset_range_max = max( offset_range );
x = ( offset_range - offset_range_mean ) / ( offset_range_max - offset_range_min);
y = normat( metric );
fitType = 'poly2';
fitobject = fit( x', y', fitType);

%% Minimum
xmin = -fitobject.p2 / 2 / fitobject.p1;
ymin = fitobject.p1*xmin^2 + fitobject.p2*xmin + fitobject.p3;
offset = ( offset_range_max - offset_range_min) * xmin + offset_range_mean;
fprintf( '\n offset minimum: %g', offset )

%% Plot metrics
figure('Name', 'OFFSET: metric');
plot( fitobject, x, y, '+');
hold on
plot( xmin, ymin ,  'o' )
axis tight
xlabel( 'offset' )
legend( {'data', 'fit', 'minimum' } )
title(sprintf('rotation axis: metric VS offset (scaled)'))
drawnow

fprintf( '\n' )