function vol = MaskingWindowSine( vol, fraction)
% Keep central volume and let boundary region decay sinusoidally to zero.
%
% vol: 2D image or 3D volume
% fraction: scalar in [0 1]. default: 0.05. fraction of the extend of the
% array over which the original array shall drop to 0.
%
% Written by Julian Moosmann, 2017-05-27. Last version:
%
% vol = MaskingWindowSine( vol, fraction)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    vol = rand( 10, 20, 8 );
end
if nargin < 2
    fraction = 0.3;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = ceil( fraction * size( vol ) );
x = s(1);
y = s(2);
z = s(3);

%% Dim 1
s = ( 1 / 2 - cos(pi * (0:x-1) / (x - 1) ) / 2)';
vol(1:x,:,:) = bsxfun(@times, vol(1:x,:,:), s );

s = ( cos( pi * (0:x-1) / (x - 1) ) / 2 + 1 / 2)';
vol(end-x+1:end,:,:) = bsxfun( @times, vol(end-x+1:end,:,:), s ); 

%% Dim 2
s = ( 1 / 2 - cos(pi * (0:y-1) / (y - 1) ) / 2);
vol(:,1:y,:) = bsxfun(@times, vol(:,1:y,:), s );

s = ( cos( pi * (0:y-1) / (y - 1) ) / 2 + 1 / 2);
vol(:,end-y+1:end,:) = bsxfun( @times, vol(:,end-y+1:end,:), s ); 

%% Dim 3
s = ( 1 / 2 - cos(pi * (0:z-1) / (z - 1) ) / 2);
vol(:,:,1:z) = bsxfun(@times, vol(:,:,1:z), shiftdim( s, -1) );

s = ( cos( pi * (0:z-1) / (z - 1) ) / 2 + 1 / 2);
vol(:,:,end-z+1:end) = bsxfun( @times, vol(:,:,end-z+1:end), shiftdim( s, -1) );
