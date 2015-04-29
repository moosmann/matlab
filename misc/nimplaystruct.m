function nimplaystruct(imstruct,field,RemoveLargeScales)
%Play sequence of images in a structure array as a video using a normalized
%input for MATLAB's implay. Optionally large scales can be removed from the
%images.
%
% imstruct: struct - struct containing the images
% field: string - name of the struct field containing the images, e.g.
%   'tie', 'tiePNLO', 'ctf', or 'ctfProjected'
% RemoveLargeScales: scalar - 0: does nothing, 1: Removes the large scales
%   form the images using default values of the 'RemoveLowFreq' function, >1:
%   uses the given input values as standard deviation for the Gaussian
%   blurring kernel in the 'RemoveLowFreq' function
%
% Written by Julian Moosmann

%% Default arguments.
if nargin<2
    field = 'tie';
end
if nargin<3
    RemoveLargeScales = 0;
end

%% Programme.
% Make an array of images out of the input struct images. !!!! Same
% variable name is used to save memory: imstuct is an array of images now.
% !!!!
eval(['imarray = cat(3,' sprintf('%s(:).%s','imstruct',field) ');']);

% Get dimension of the image array.
[dimX, dimY, dimZ] = size(imarray);

% Remove the image's large scales by subtraction.
if RemoveLargeScales == 1
    for ii = dimZ:-1:1
        imarray(:,:,ii) = RemoveLowFreq(imarray(:,:,ii));
    end
elseif RemoveLargeScales > 1
    for ii = dimZ:-1:1
        imarray(:,:,ii) = RemoveLowFreq(imarray(:,:,ii),RemoveLargeScales);
    end
end

% Find minimum and maximum of each matrix in the input array and create a
% an array corresponding the input array to subtract the values from the
% input arrray.
imarrayMin = repmat(min(min(imarray)),[dimX,dimY,1]);
imarrayMax = repmat(max(max(imarray)),[dimX,dimY,1]);

% Renormalize: subtract minimum, then divide by maximum-minimum.
imarray    = (imarray-imarrayMin)./(imarrayMax-imarrayMin);

% Pay images using implay.
implay(imarray)