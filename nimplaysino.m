function nimplaysino(imstack,normGlobal)
%Play stack of images (array) as video clip. Contrast is adapted by
%normalizing the images before. Images have to be stacked along the third
%dimension.  
%
%nimplay(imstack,normGlobal)

%% Default arguments
if nargin < 2
    normGlobal = 1;
end

imstack = permute(imstack,[3 2 1]);
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
implay(imstack);
