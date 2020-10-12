function FTofIm = ftool(im,epsilon,subtractMean,doFFTshift,PaddingFactor, fftdir)
%Show logarithm of modulus of Fourier-transformed image 'im'. Mean is
%substracted before processing. If 'im' is an image stack, 'slice' picks
%the correspondig image from the stack. 'espilon' removes a possible
%singularity at zero frequency when the logarithm is taken:
%log(10^(-epsilon)+abs(fftshift(fft2(im)))) 
%
% ftool(im,epsilon,subtractMean,doFFTshift,PaddingFactor)
%
% Written by Julian Moosmann, last version: 2013-10-24

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    epsilon = 0;
end
if nargin < 3
    subtractMean = false;
end
if nargin < 4
    doFFTshift = true;
end
if nargin < 5
    PaddingFactor = 0;
end
if nargin < 6
    fftdir = 'both';
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = double(squeeze(im));
%% Image parameters
imMin  = min(im(:));
imMax  = max(im(:));
imMean = mean(im(:));
imStd  = std(im(:));
[dimx,dimy] = size(im);
fprintf('Size: %u x %u. ',dimx,dimy)
%% Padding
if PaddingFactor > 0
    im = padarray( im, [ PaddingFactor * dimx, PaddingFactor * dimy ], 'symmetric', 'post');    
    fprintf('\nPadded size: %u x %u. ',size(im)) 
end
fprintf( '\nepsilon: %f', epsilon )
%% Subtract mean
if subtractMean > 0
    fprintf('\nSubstract mean: %g',imMean)
    im = im - mean(im(:));
else
    fprintf('\nNo substraction of mean')
end
%% FT of image
switch fftdir
    case 'both'
        im = fft2(im);
    case {'first', '1st', 1}
        im = fft(im, [], 1);
    case{'second', '2nd', 2}
        im = fft(im, [], 2);
end
%% fftshift
if doFFTshift > 0
    fprintf('\nApply fftshift')
    im = fftshift(im);
    fprintf('\nShow: log( 10^(-%g) + abs( fftshift( fft2(im) ) ) )',epsilon)
else
    fprintf('\nNo fftshift')
    fprintf('\nShow: log( 10^(-%g) + abs( fft2(im) ) )',epsilon)
end
if ~PaddingFactor
    fprintf('\nNo Padding. ')
end
%% Output
if nargout > 0
    FTofIm = im;
end
%% Regularized logarithim of modulus of FT of image
im = log(10^(-epsilon)+abs(im));
%% Print info
fprintf('Input image: [Min Max Min-Max Mean Std] = [%g %g %g %g %g]\n',imMin,imMax,imMax-imMin,imMean,imStd);
imMin  = min(im(:));
imMax  = max(im(:));
imMean = mean(im(:));
imStd  = std(im(:));
fprintf('FT of image: [Min Max Min-Max Mean Std] = [%g %g %g %g %g]\n',imMin,imMax,imMax-imMin,imMean,imStd);
%% Show FT of image
itool(im,inputname(1))
