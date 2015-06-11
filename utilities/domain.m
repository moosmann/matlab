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
fprintf(' %s: Min=%g, Max=%g, Mean=%g, Sta=%g, Var=%g Max-Min=%g\n', ...
    imNameString,imMin,imMax,imMean,imSta,imVar,imMax-imMin);
