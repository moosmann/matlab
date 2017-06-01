function arr = MaskingExponential( arr, fraction)
% Keep central ellipse within the 'radial_fraction' of the largest ellipse
% fitting within the rectangular image and set the region outside of this
% ellipse to the mean within the disc or to 'value'.
%
% arr: 2D image or 3D volume
% fraction: scalar in [0 1]. default: 0.05. fraction of the extend of the
% array over which the original array shall drop to 0.
%
% Written by Julian Moosmann, 2017-05-27. Last version:
%
% arr = MaskingExponential( arr, fraction)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    arr = rand( 10, 20, 8 );
end
if nargin < 2
    fraction = 0.3;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = ceil( fraction * size( arr ) );
x = s(1);
y = s(2);
z = s(3);

%% Dim 1
s = ( 1 / 2 - cos(pi * (0:x-1) / (x - 1) ) / 2)';
arr(1:x,:,:) = bsxfun(@times, arr(1:x,:,:), s );

s = ( cos( pi * (0:x-1) / (x - 1) ) / 2 + 1 / 2)';
arr(end-x+1:end,:,:) = bsxfun( @times, arr(end-x+1:end,:,:), s ); 

%% Dim 2
s = ( 1 / 2 - cos(pi * (0:y-1) / (y - 1) ) / 2);
arr(:,1:y,:) = bsxfun(@times, arr(:,1:y,:), s );

s = ( cos( pi * (0:y-1) / (y - 1) ) / 2 + 1 / 2);
arr(:,end-y+1:end,:) = bsxfun( @times, arr(:,end-y+1:end,:), s ); 

%% Dim 3
s = ( 1 / 2 - cos(pi * (0:z-1) / (z - 1) ) / 2);
arr(:,:,1:z) = bsxfun(@times, arr(:,:,1:z), shiftdim( s, -1) );

s = ( cos( pi * (0:z-1) / (z - 1) ) / 2 + 1 / 2);
arr(:,:,end-z+1:end) = bsxfun( @times, arr(:,:,end-z+1:end), shiftdim( s, -1) );



