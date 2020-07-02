function imf = FilterLine(im,theta,width,ShowInfo)

%% Default arguments
if nargin < 2
    theta = 76.15;% angle from the horizontal in degrees
end
if nargin < 3
    width = 4;% width in pixel
end
if nargin < 4
    ShowInfo = 1;
end

%% Body
%im = (im0(1:2:end,1:2:end)+im0(2:2:end,1:2:end)+im0(1:2:end,2:2:end)+im0(2:2:end,2:2:end))/4;
[dimx ,dimy] = size(im);
% Create mesh
[x ,y] = meshgrid(-dimy/2:dimy/2-1,-dimx/2:dimx/2-1);
% Create filter
m=ones(dimx,dimy);
m(round(abs(x-cotd(theta)*y))<width) = 0;
m=fftshift(m);
% Filter image
imf = real(ifft2(m.*fft2(im)));
%% Show info
if ShowInfo == 1
    ishow(m)
    ishow(im-imf)
    ishow(im)
    ishow(imf)
end
