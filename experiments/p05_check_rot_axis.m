function p05_check_rot_axis(im1, im2, rotate_images )
% Compute tilt of images that are projections which differ by a rotation
% angle of pi.

%% Default arguments
if nargin < 1
    % Read images
    par_folder = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/images/';
    rows = [1200 2200];
    cols = [1 3056];    
    im1 = FilterPixel( flipud( imread([par_folder 'holder_000'],'tif', 'PixelRegion', {rows, cols}) ), [0.01, 0.001] );
    im2 = FilterPixel( flipud( imread([par_folder 'holder_180'],'tif', 'PixelRegion', {rows, cols})), [0.01, 0.001] );
    % Normalize
    im1 = im1 / mean( im1(:) );
    im2 = im2 / mean( im2(:) );
end
if nargin < 3
    rotate_images = 0;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rotate_images
    im1 = rot90(im1);
    im2 = rot90(im2);
end


im0 = im1;
im180 = im2;


%h = figure('Name', 'proj 0 degree');
subplot(2,2,1)
imsc( im0)
axis equal tight
title('0 degrees')
subplot(2,2,2)
imsc( im180)
title('180 degrees')
axis equal tight

out = ImageCorrelation(im0, fliplr(im180), 0);
rot_axis_offset =  out.Yshift / 2 ;
rot_axis_pos = out.VerticalRotationAxisPosition;
fprintf('\n rotation axis position: %g', rot_axis_pos)
fprintf('\n rotation axis offset: %g', rot_axis_offset)

im0c = RotAxisSymmetricCropping( im0, rot_axis_pos, 2);
im180c = fliplr(RotAxisSymmetricCropping( im180, rot_axis_pos, 2) );

tilt_max = 0.5;
tilt_stride = 0.05;
tilt_range = 0:tilt_stride:tilt_max;
tilt_range = [-tilt_range(end:-1:2), tilt_range];

err_l2 = zeros( 1, numel( tilt_range) );
err = zeros( [size(im0c), numel( tilt_range) ]);
for nn = 1:numel(tilt_range)
    tilt = tilt_range(nn);
    im0r = imrotate(im0c, -tilt/2, 'bilinear', 'crop');
    im180r = imrotate(im180c, tilt/2, 'bilinear', 'crop');
    d = abs(im0r - im180r);
    err(:, :, nn) = d;
    err_l2(nn) = sqrt( mean2( d(201:end-200,201:end-200) ).^2 ) ;
end

subplot(2,2,3)
[~, tilt_opt_pos] = min( err_l2 );
tilt_opt = tilt_range(tilt_opt_pos);
imsc( err(:, :, tilt_opt_pos) )
title('error map at optimal tilt')
fprintf('\n optimal tilt angle: %g\n', tilt_opt ) 

subplot(2,2,4)
plot(tilt_range, err_l2)
axis tight
title( 'L2 norm of error map' )

im0r = imrotate(im0c, -tilt_opt/2, 'bilinear', 'crop');
im180r = imrotate(im180c, tilt_opt/2, 'bilinear', 'crop');

outr = ImageCorrelation(im0r, im180r,0);
rot_axis_offsetr =  outr.Yshift ;
rot_axis_posr = outr.VerticalRotationAxisPosition;
fprintf('\n rotation axis position after tilt correction: %g', rot_axis_posr)
fprintf('\n rotation axis offset after tilt correction: %g', rot_axis_offsetr)



%nimplay( cat(3, im0c, im180c ) )

fprintf('\n')