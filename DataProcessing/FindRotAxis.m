function [offset, metric] = FindRotAxis( sino, shift, sampleWidth_to_detectorWidth )
% Determine the center of rotation via Fourier analysis of projection data.
% See articles:
%
% Written by J.Moosmann,2018-08-05. Last version: 2018-08-05

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    phan = phantom( 128 );
    sino = radon( phan );
end
if nargin < 2
    shift = 0;
end
if nargin < 3
    sampleWidth_to_detectorWidth = 1;
end
if nargin < 4
    more_figs = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove lateral sample shift
sino = CropShift( sino, shift );
[s1, s2] = size( sino );
% Offset search range
offset_range = (1:round( 1.5 * s1 / 10 ):s1) - floor( s1 / 2 ) ;

% Triangular mask
% Fourier coordinates 1D
precision = 'double';
xi  = FrequencyVector(2*s2,precision,1);
eta = FrequencyVector(2*s1,precision,1);
% Fourier coordinates 2D
[xi, eta]  = meshgrid(xi,eta);
mask = eta - abs( sampleWidth_to_detectorWidth * xi ) < 0;
mask = mask & flipud( mask );
sum_mask = sum( mask(:) );
% Loop over offset range
metric = zeros( size( offset_range ) );
if more_figs
    h0 = figure( 'Name', 'mask' );
    figure( h0)
    imsc( fftshift( mask ) )
    drawnow
    h1 = figure( 'Name', 'sino and mask .* FT sino' );
end

for nn = 1:numel( offset_range )
    offset = offset_range(nn);
    if offset > 0
        tmp = cat( 2, sino, cat( 1, sino(end-offset+1:end,:), sino(1:end-offset,:)) );
    else
        tmp = cat( 2, sino, cat( 1, sino(1-offset:end,:), sino(1:-offset,:)) );
    end
    if more_figs
        figure( h1)
        subplot(1,2,1)
        imsc( tmp )
        title('sinogram')%drawnow;pause(1);
        axis tight equal
    end
    tmp = padarray( tmp, [s1, 0], 'symmetric', 'post' );
    tmp = fft2( tmp );
    tmp = abs( tmp );
    tmp = mask .* tmp;
    if more_figs
        subplot(1,2,2)
        imsc( fftshift( log(1+abs(tmp) ) ) )
        title('mask .* FT sino')
        axis tight equal
        drawnow;pause(1);
    end
    metric(nn) = norm( tmp, 1 ) / sum_mask;
end

metric = normat( metric );
offset = min( metric(:) );

%% Fit
offset_range_mean = mean( offset_range );
offset_range_min = min( offset_range );
offset_range_max = max( offset_range );
x = ( offset_range - offset_range_mean ) / ( offset_range_max - offset_range_min);
y = metric;
fitType = 'poly2';
fitobject = fit( x', y', fitType);

%% Plot
% Plot metrics
figure('Name', 'OFFSET: metric');
plot( fitobject, x, metric, '+');
axis tight
xlabel( 'offset' )
%legend( metrics_offset(x).name )
%ax1 = gca;
%ax2 = axes( 'Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
%line(1:numel( offset ), 0, 'Parent', ax2 )
%xlabel( 'index (image no.)' )
%set( ax1, 'YTick', [] )
%set( ax2, 'YTick', [] )
title(sprintf('rotation axis: metrics VS offset'))
drawnow

