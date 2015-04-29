function out = fshow(im,slice,epsilon,padding)
%Show logarithm of modulus of Fourier-transformed image 'im'. Mean is
%substracted befor processing. If 'im' is an image stack, slice pick the
%correspondig image out of the stack. 'espilon' removes possible
%singularity when logarithm is taken: log(10^(-epsilon)+abs(fftshift(fft2(im))))

%% Default arguments
if nargin < 2
    slice = 1;
end;
if nargin < 3
    epsilon = 0;
end
if nargin < 4
    padding = 0;
end
%% Pick slice if im is a stack, subtract mean, FT image, take modulus and logarithm.
im = squeeze(im);
im = im(:,:,slice);
im = im - mean(im(:));
if padding > 0
    if mod(size(im,1),2)
        im = im(1:end-1,:);
    end
    if mod(size(im,2),2)
        im = im(:,1:end-1);
    end
    im = padarray(im,size(im)/2,'symmetric','both');
    im = fftshift(fft2(im));
    im = im(1:2:end,1:2:end)+im(2:2:end,1:2:end)+im(1:2:end,2:2:end)+im(2:2:end,2:2:end);
    im = log10(10^(-epsilon)+abs(im));
else
    im = log10(10^(-epsilon)+abs(fftshift(fft2(im))));
end
%% Show logarithm of modulus of FTed image.
domain(im,1,sprintf('log(%g+abs(fftshift(fft2(im))))',10^epsilon))
figure
imshow(im,[],'InitialMagnification','fit')
colorbar;
%% output
if nargout == 1
    out = im;
end