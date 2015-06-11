function imageArray = norar(imageArray)
%Renormalize 3D input imageArrayray 'imageArray' to the interval [0,1].

[d1,d2,~] = size(imageArray);

% Find minimum and maximum of each matrix in the input imageArrayray and create a
% an imageArrayray corresponding the input imageArrayray to subtract the values from the
% input imageArrayrray.
imageArraymin = repmat(min(min(imageArray)),[d1,d2,1]);
imageArraymax = repmat(max(max(imageArray)),[d1,d2,1]);

% Renormalize: subtract minimum, then divide by maximum-minimum.
imageArray     = (imageArray-imageArraymin)./(imageArraymax-imageArraymin);