
function FTofIm = ftool(im,epsilon,subtractMean,doFFTshift,PaddingFactor)
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
    im = padarray(im,[ ceil(PaddingFactor*dimx/2), ceil(PaddingFactor*dimy/2)],'symmetric','pre');
    im = padarray(im,[floor(PaddingFactor*dimx/2),floor(PaddingFactor*dimy/2)],'symmetric','post');
    fprintf('Padded size: %u x %u. ',size(im)) 
end
%% Subtract mean
if subtractMean > 0
    fprintf('Substract mean: %g. ',imMean)
    im = im - mean(im(:));
else
    fprintf('No substraction of mean. ')
end
%% FT of image
im = fft2(im);
%% fftshift
if doFFTshift > 0
    fprintf('Do fftshift.\n')
    im = fftshift(im);
    fprintf('Show: log( 10^(-%g) + abs( fftshift( fft2(im) ) ) )\n',epsilon)
else
    fprintf('No fftshift.\n')
    fprintf('Show: log( 10^(-%g) + abs( fft2(im) ) )\n',epsilon)
end
if ~PaddingFactor
    fprintf('No Padding. ')
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
