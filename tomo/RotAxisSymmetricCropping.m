function sino = RotAxisSymmetricCropping(sino,rotationAxisPosition)
% Crop sino symmetrically around rotation axis, such that cropped image has
% an even number of pixels and that the position of the rotation axis is
% at size(sino,2)/2.
%
% Written by Julian Moosmann, last version 2013-10-21

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rot axis is centered
if rotationAxisPosition == size(sino,2)/2
    return
% Rot axis is left to center    
elseif rotationAxisPosition < size(sino,2)/2
    sino = sino(:,1:round(2*rotationAxisPosition),:);
% Rot axis is right to center
else
    horWidth = size(sino,2)-rotationAxisPosition;
    sino = sino(:,end-round(2*horWidth)+1:end,:);
end