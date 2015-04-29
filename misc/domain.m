function domain(im,slice,imNameString)
% Print range [minimum maximum], mean, standard deviation, variance, maximum-minimum.

%% Default arguments.
if nargin < 2
    slice = 1;
end
if nargin < 3
    imNameString = inputname(1);
end;
%% Main code
im = im(:,:,slice);
if  ~isa(im,'single')
    im     = double(im);
end
imMin  = min(real(im(:)));
imMax  = max(real(im(:)));
imMean = mean(im(:));
imVar  = var(im(:));
imSta  = sqrt(imVar);
%% Print info
fprintf(' %s: Min=%8g, Max=%8g, Mean=%8g, Sta=%8g, Var=%8g Max-Min=%8g\n', ...
    imNameString,imMin,imMax,imMean,imSta,imVar,imMax-imMin);
