function domains(varargin)
% Print range [minimum maximum], mean, standard deviation, variance, maximum-minimum.

%%
slice = 1;
%% Main code
% Loop over input variables
for nn=1:numel(varargin)
    im = varargin{nn};
    im = im(:,:,slice);
    im     = double(im);
    imMin  = min(real(im(:)));
    imMax  = max(real(im(:)));
    imMean = mean(im(:));
    imVar  = var(im(:));
    imSta  = sqrt(imVar);
    %% Print info
    fprintf(' %s: Min=%g, Max=%g, Mean=%g, Sta=%g, Var=%g Max-Min=%g\n', ...
        inputname(nn),imMin,imMax,imMean,imSta,imVar,imMax-imMin);
end