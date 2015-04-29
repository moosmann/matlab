function [imf DiskFilter ] = FilterDisk(im,DiskFilterWidth,RetainMean,DoFFTshift,DoShowInfo,Padding)
% Cross and center-symmetric filter with Gaussian profile for high-pass and
% directional filtering along x- and y-axes of images in Fourier space

%% Default arguments
if nargin < 2
    DiskFilterWidth = 60;
end
if nargin < 3
    RetainMean = 1;
end
if nargin < 4
    DoFFTshift = 1;
end
if nargin < 5
    DoShowInfo = 1;
end
if nargin < 6
    Padding = 1;
end
%% Body
%im = (im0(1:2:end,1:2:end)+im0(2:2:end,1:2:end)+im0(1:2:end,2:2:end)+im0(2:2:end,2:2:end))/4;
%% Padding
[dim1 dim2] = size(im);
if Padding > 1
    im = padarray(im,[dim1/2 dim2/2],'symmetric','both');
end
[dimx dimy] = size(im);
% Create mesh
[x y] = meshgrid(-dimy/2:dimy/2-1,-dimx/2:dimx/2-1);
%% Circular center-symmetric high-pass filter in Fourier space
DiskFilter = 1-exp(-(x.^2+y.^2)/DiskFilterWidth^2/2);
if RetainMean == 1
    DiskFilter(dimx/2+1,dimy/2+1) = 1;% Restore mean
end
if DoFFTshift == 1
    DiskFilter = fftshift(DiskFilter);
end
%% Apply filter to image
if DoFFTshift == 1
    imf = real(ifft2(DiskFilter.*fft2(im)));
else
    imf = real(ifft2(fftshift(DiskFilter.*fftshift(fft2(im)))));
end
%% Unpadding
if Padding > 1
    imf = imf(dim1/2+(1:dim1),dim2/2+(1:dim2));
end
%% Show info
if DoShowInfo == 1
    domain(DiskFilter);
    domain(im);domain(imf)
    if DoFFTshift == 0
        ishow(DiskFilter);
    else
        ishow(fftshift(DiskFilter));
    end
    ishow(im);
    ishow(imf)
end


