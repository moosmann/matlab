function nimplay(imstack,normGlobal,permuteOrder, name)
% Play stack of images (array) as video clip. Contrast is adapted by
% normalizing the images before hand. Images have to be stacked along the
% third dimension, otherwise use permuteOrder to rearrange stack dimensions.
% 
% imstack: 3D-array of images: vert x hor x slice (Matlab notation)
% normGlobale: scalar, default: 0. 1: subtract from each slice its mean. 0:
% subtract from each slice the gobal mean.
% permuteOrder: vector. Permute array dimension.
%
% Written by Julian Moosmann, last version 2013-11-01
%
%nimplay(imstack,normGlobal,permuteOrder)

%% Default arguments
if nargin < 2
    normGlobal = 0;
end
if nargin < 3
    permuteOrder = 0;
end
if nargin < 4
    name = '';
end

%% Permute stack
if permuteOrder(1)
    imstack = permute(imstack,permuteOrder);
end
%% Compute minima and maxima
if normGlobal
    armin = min(imstack(:));
    armax = max(imstack(:));
else
    % Find minimum and maximum of each matrix in the input array and create a
    % an array corresponding the input array to subtract the values from the
    % input arrray.
    [d1, d2, ~] = size(imstack);
    armin = repmat(min(min(imstack)),[d1,d2,1]);
    armax = repmat(max(max(imstack)),[d1,d2,1]);
end
%% Renormalize: subtract minimum, then divide by maximum-minimum.
imstack = (imstack-armin)./(armax-armin);
%% Play normalized image array as a clip.
h = implay(imstack);
if isempty( name )
    name = sprintf( 'volume shape: %u x %u x %u', size(imstack) );
end
set(h.Parent, 'Name', name)
