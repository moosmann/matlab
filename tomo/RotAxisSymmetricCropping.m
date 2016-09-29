function sino = RotAxisSymmetricCropping(sino,rotationAxisPosition, direction)
% Crop sino symmetrically around rotation axis, such that cropped image has
% an even number of pixels and that the position of the rotation axis is
% at size(sino,2)/2.
%
% Written by Julian Moosmann, last version 2013-10-21

if nargin < 3
    direction = 2;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rot axis is centered
if rotationAxisPosition == size( sino, direction ) / 2
    return
% Rot axis is left to center    
elseif rotationAxisPosition < size( sino, direction ) / 2
    if direction == 1
        sino = sino(1:round(2*rotationAxisPosition), :, :);
    elseif direction == 2
        sino = sino(:, 1:round(2*rotationAxisPosition), :);
    end
% Rot axis is right to center
else
    newWidth = size( sino, direction) - rotationAxisPosition;
    if direction == 1
        sino = sino(end-round(2*newWidth)+1:end, :, :);        
    elseif direction == 2
        sino = sino( :, end-round(2*newWidth)+1:end, :);
    end
end