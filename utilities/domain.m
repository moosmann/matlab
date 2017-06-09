function domain(im,slice,imNameString)
% Print range [minimum maximum], mean, standard deviation, variance, maximum-minimum.

%% Default arguments.
if nargin >= 2
    im = im(:,:,slice);    
end
if nargin < 3
    imNameString = inputname(1);
end;
%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  ~isa(im,'single')
    im     = double(im);
end
imMin  = min(real(im(:)));
imMax  = max(real(im(:)));
imMean = mean(im(:));
imVar  = var(im(:));
imSta  = sqrt(imVar);
%% Print info
fprintf(' %s: Min=%12g, Max=%12g, Mean=%12g, Sta=%12g, Var=%12g Max-Min=%12g\n', ...
    imNameString,imMin,imMax,imMean,imSta,imVar,imMax-imMin);
