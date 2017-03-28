function [out, corr_map] = ImageCorrelation2(im1, im2, padding, use_gpu, verbose)
%Find the relative movement of image im2 w.r.t. to im1 by means of
%correlating the images. Thus rotation axis can be determined. Allows for
%subpixel precsision.
%
% Written Julian Moosamnn. Last modification: 2017-02-20, GPU support added

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    im1 = zeros(8);
    im1(4,3:5) = [1 3 2];
end
if nargin < 2
    im2 = zeros(8);
    im2(4,4:6) = [2 5 1];
end
if nargin < 3
    padding = 0;
end
if nargin < 4
    use_gpu = 0;
end
if nargin < 5
    verbose = 1;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Padding
im1 = padarray( im1, padding * size( im1 ), 0, 'post');
im2 = padarray( im2, padding * size( im2 ), 0, 'post');

% Normalize
im1 = SubtractMean( im1 ) / std2( im1 );
im2 = SubtractMean( im2 ) / std2( im2 ) ;

% Image shape and center
[d1, d2] = size(im1);
cent1 = (d1 + 1) / 2 ;
cent2 = (d2 + 1) / 2 ;

% GPU
if use_gpu
    im1 = gpuArray( im1 );
    im2 = gpuArray( im2 );
end

res = 1;

% Cross-correlation map
corr_map = fft2( im1 ) .* fft2( rot90(im2, 2) );
corr_map = padarray( corr_map, (res - 1) * size(corr_map), 0, 'post');
corr_map = fftshift( ifft2( corr_map, 'symmetric' ) );
%whos corr_map

% Correlation map maximum
val = max( corr_map(:) );
[m2, m1] = meshgrid( 1:d2, 1:d1);
corr_map = corr_map == val;
N = sum( corr_map(:));
pos1 = corr_map .* m1;
pos1 = sum( pos1(:) ) / N;
pos2 = corr_map .* m2;
pos2 = sum( pos2(:) ) / N;
fprintf( '\n pos1 : %g \n pos2 : %g \n', pos1, pos2 )

%fprintf( '\n\n pos2 : %g \n\n', pos2)
pos1 = - pos1 + res * d1 / 2;
pos2 = - pos2 + res * d2 / 2;

%% Check
% Add 0.5 or 1 pixel to peak position for even or odd dimensions
% switch mod(d1,2)
%     case 0
%         offsetX = 0.5;
%     case 1
%         offsetX = 1;
% end
% switch mod(d2,2)
%     case 0
%         offsetY = 0.5;
%     case 1
%         offsetY = 1;
% end
% pos1 = pos1 + offsetX;
% pos2 = pos2 + offsetY;

% Relative image shift
out.shift1 = gather( pos1 );
out.shift2 = gather( pos2 );

% Rotation axis position
out.rot_axis_pos_1 = gather( cent1 + pos1 );
out.rot_axis_pos_2 = gather( cent2 + pos2 );

imsc( corr_map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print results
if verbose
    fprintf('image shape : %u x %u\n', d1, d2);
    fprintf('image center : [%g %g]\n', cent1, cent2);
    fprintf('correlation map peak : [%.2f %.2f], value: %g\n', pos1, pos2, val);
    fprintf('relative shift : [%f %f]\n', out.shift1, out.shift2);
    fprintf('rotation axis pos : [%f %f]\n', out.rot_axis_pos_2, out.rot_axis_pos_1);
end
