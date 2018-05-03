function nimplay( vol, renorm_slicewise, permute_order, figure_name)
% Play 3D array using MATLAB's movie play. Contrast is adjuste slicewise by
% default unless renorm_slicewise is set to zero. Movie is run along 3rd
% dimension unless dimension are rearranged using the permutation order
% parameter.
% 
% vol : 3D-array of images
% renorm_slicewise : scalar, default: 1. Adjust colorbar slicewise, losing scaling
%  information relative to other slices. If > 0 and < 1: adjust to
%  renorm_slicewise * min/max.
% permute_order : empty or 3D vector. Permute array dimension. If empty do
%  not permute.
% figure_name : string (of figure)
%
% Written by Julian Moosmann, last version 2013-11-01
% Last modification: 2017-10-09
%
%nimplay( vol, renorm_slicewise, permute_order, figure_name)

%% Default arguments
if nargin < 2
    renorm_slicewise = 1;
end
if nargin < 3
    permute_order = [];
end
if nargin < 4
    figure_name = '';
end

%% Permute dimensions
if ~isempty( permute_order ) && prod( permute_order ) ~= 0
    vol = permute(vol, permute_order);
end

if isinteger( vol )
    vol = single( vol );
end

%% Normalization
if renorm_slicewise == 0
    armin = min( vol(:) );
    armax = max( vol(:) );
else
    % Find minimum and maximum of each matrix in the input array and create a
    % an array corresponding the input array to subtract the values from the
    % input arrray.
    [d1, d2, ~] = size(vol);
    armin = renorm_slicewise * repmat(min(min(vol)),[d1,d2,1]);
    armax = renorm_slicewise * repmat(max(max(vol)),[d1,d2,1]);
end

%% Renormalize: subtract minimum, then divide by maximum-minimum.
vol = ( double( vol ) - armin ) ./ ( armax - armin );

%% Play normalized volume in movie player
h = implay(vol);
if isempty( figure_name )
    figure_name = sprintf( 'volume shape: %u x %u x %u', size(vol) );
end
set(h.Parent, 'name', figure_name)
drawnow
pause(0.01)
