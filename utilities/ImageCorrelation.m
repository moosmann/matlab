function [out, map] = ImageCorrelation(im1, im2, verbose, show_plots, use_gpu, use_gradient, blur_images)
%Find the relative movement of image im2 w.r.t. to im1 by means of
%correlating the images. Thus rotation axis can be determined. Allows for
%subpixel precsision using the center of mass of a samll region centered
%around the maximum of the correlation map.
%
% Written Julian Moosamnn. Last modification: 2017-02-22, GPU support added

% TODO: fix plotting option if axis is close to boundary and center of mass
% region is too large

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    nn = 21;
    im1 = zeros(nn);
    im1(10,9) = 1;
end
if nargin < 2
    im2 = zeros(nn);
    im2(10,10:11) = 1;
end
if nargin < 3
    verbose = 0;
end
if nargin < 4
    show_plots = 0;
end
if nargin < 5
    use_gpu = 0;
end
if nargin < 6
    use_gradient = 0;
end
if nargin < 7
    blur_images = 0;
end
%% Notes
% Normalize: Don't do it. Make sure to have physical values in your
% intensity map. Normalization does mess up center of mass computation.

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum radius of region centered around the correlation map peak where
% the centere of mass is calculated on
radius = 5;

% Image center
[d1, d2] = size(im1);
cen_1 = ( d1 + 1 ) / 2;
cen_2 = ( d2 + 1 ) / 2;

% GPU
if use_gpu
    im1 = gpuArray( im1 );
    im2 = gpuArray( im2 );
end

% Blur
if blur_images
    im1 = FilterBlur( im1, [3 3], 10);
    im2 = FilterBlur( im2, [3 3], 10);
end

%
if use_gradient
    l = 2;
    [g1, g2] = gradient(im1);
    im1 = abs(g1).^l + abs(g2).^l;
    [g1, g2] = gradient(im2);
    im2 = abs(g1).^l + abs(g2).^l;
end

% Correlation map
map = fft2( im1 ) .* fft2( rot90(im2,2) );
map = ifft2( map, 'symmetric' );
map = fftshift( map );
map = abs( map );

% Find the value and the (index) position of the maximum of the correlation map.
[val, ind] = max( map(:) );

% Convert linear index to row and column subscripts
[pos_1, pos_2] = ind2sub( [d1 d2], ind);

% Add 0.5 or 1 pixel to peak position for even or odd dimensions
switch mod( d1, 2)
    case 0
        offset_1 = 0.5;
    case 1
        offset_1 = 1;
end
switch mod( d2, 2)
    case 0
        offset_2 = 0.5;
    case 1
        offset_2 = 1;
end

% Minimum radius
radius_1 = min( [radius, pos_1 - 1, d1 - pos_1]);
radius_2 = min( [radius, pos_2 - 1, d2 - pos_2]);

% ROI centered at correlation map peadk for COM computation
roi = map(pos_1 + (-radius_1:radius_1), pos_2 + (-radius_2:radius_2));

% x-/y-coordinate-matrices 
[mdy, mdx]   = meshgrid(1:d2,1:d1);

% Total mass
M = sum(roi(:));

% Center of mass of ROI
com_1 = sum(sum( mdx(pos_1 + (-radius_1:radius_1), pos_2 + (-radius_2:radius_2)) .* roi)) / M + offset_1;
com_2 = sum(sum( mdy(pos_1 + (-radius_1:radius_1), pos_2 + (-radius_2:radius_2)) .* roi)) / M + offset_2;

%pos_1 = pos_1 + offset_1;
%pos_2 = pos_2 + offset_2;

% Relative shift
out.shift1 = gather( com_1 - cen_1 );
out.shift2 = gather( com_2 - cen_2 );

% Rotation axis position
out.rot_axis_pos1 = gather( cen_1/2 + com_1/2 );
out.rot_axis_pos2 = gather( cen_2/2 + com_2/2 );

%% Print  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    fprintf('image shape : %u x %u\n', d1, d2);
    fprintf('image center : [%g %g].\n', cen_1, cen_2);
    fprintf('correlation map peak : [%.2f %.2f], value: %g\n', pos_1, pos_2, val);
    fprintf('center of mass of ROI within radi [%u %u] centered at peak: [%f %f]\n', radius_1, radius_2, com_1, com_2);
    fprintf('relative shift : [%f %f].\n',out.shift1, out.shift2);
    fprintf('rotation axis position : [%f, %f]\n', out.rot_axis_pos1, out.rot_axis_pos2);
end
%% Show plots 
if show_plots    
    map = gather( abs( map ) );    
    h = imtool( map, [], 'InitialMagnification', 'fit');
    set( h, 'Name', sprintf( 'Correlation map. peak : %u, %u. COM radi : %u, %u', pos_1, pos_2, radius_1, radius_2 ) )
end
