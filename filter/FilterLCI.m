function [imf CircularFilter CrossFilter] = FilterLCI(im,CircularFilterWidth,CrossFilterWidth,RetainMean,DoFFTshift,DoShowInfo)
% Cross and center-symmetric filter with Gaussian profile for high-pass and
% directional filtering along x- and y-axes of images in Fourier space

%% Default arguments
if nargin < 2
    CircularFilterWidth = 60;
end
if nargin < 3
    CrossFilterWidth = 2;
end
if nargin < 4
    RetainMean = 1;
end
if nargin < 5
    DoFFTshift = 1;
end
if nargin < 6
    DoShowInfo = 1;
end

%% Body
%im = (im0(1:2:end,1:2:end)+im0(2:2:end,1:2:end)+im0(1:2:end,2:2:end)+im0(2:2:end,2:2:end))/4;
[dimx dimy] = size(im);
% Create mesh
[x y] = meshgrid(-dimy/2:dimy/2-1,-dimx/2:dimx/2-1);
%% Circular center-symmetric high-pass filter in Fourier space
CircularFilter = 1-exp(-(x.^2+y.^2)/CircularFilterWidth^2/2);
if RetainMean == 1
    CircularFilter(dimx/2+1,dimy/2+1) = 1;% Restore mean
end
if DoFFTshift == 1
    CircularFilter = fftshift(CircularFilter);
end
%% Directional cross-like filter in Fourier space
% w1x = 1- window(@gausswin,dimx,CrossFilterWidth);
% w1y = 1- window(@gausswin,dimy,CrossFilterWidth);
% [wx wy] = meshgrid(w1y,w1x);
% CrossFilter = (wx.*wy);
CrossFilter = (1-exp(-x.^2/CrossFilterWidth^2/2)).*(1-exp(-y.^2/CrossFilterWidth^2/2));
if RetainMean == 1
    CrossFilter(dimx/2+1,dimy/2+1) = 1;% Restore mean
end
if DoFFTshift == 1
    CrossFilter = fftshift(CrossFilter);
end
%% Apply filter to image
if DoFFTshift == 1
    imf = real(ifft2(CrossFilter.*CircularFilter.*fft2(im)));
else
    imf = real(ifft2(fftshift(CrossFilter.*CircularFilter.*fftshift(fft2(im)))));
end
%% Show info
if DoShowInfo == 1
    domain(CircularFilter);
    domain(CrossFilter);
    domain(im);domain(imf)
    if DoFFTshift == 0
        ishow(CircularFilter);
        ishow(CrossFilter);
        ishow(CircularFilter.*CrossFilter);
    else
        ishow(fftshift(CircularFilter));
        ishow(fftshift(CrossFilter));
        ishow(fftshift(CircularFilter.*CrossFilter));
    end
    ishow(im);ishow(imf)
end


