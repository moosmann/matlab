function [imMin imMax imMean imSta imVar imMaxMin] = Domain(im,slice,imNameString,printInfo)
% Print range [minimum maximum], mean, standard deviation, variance,
% maximum-minimum, and output these values.
%% Default arguments.
if nargin < 2
    slice = 1;
end
if nargin < 3
    imNameString = inputname(1);
end;
if nargin < 4
    printInfo = 1;
end
%% Main code
im = im(:,:,slice);
im     = double(im);
imMin  = min(real(im(:)));
imMax  = max(real(im(:)));
imMean = mean(im(:));
imVar  = var(im(:));
imSta  = sqrt(imVar);
imMaxMin = imMax-imMin;
%% Print info
if printInfo == 1
    fprintf('Domain of %s: [%g,%g], Mean=%g, Sta=%g, Var=%g Max-Min=%g\n', ...
        imNameString,imMin,imMax,imMean,imSta,imVar,imMaxMin);
end
