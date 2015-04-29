

DataPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/test/data/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/';
HPT = 0.04;
%% Read data image
im   = FilterHotPixel(double(imread([DataPath 'proj_00551.tif'])),HPT);
%% Create Fourier space filter
% Gaussian filter
GW = 150;
[dimx dimy] = size(im);
w1x = 1- window(@gausswin,dimx,GW);
w1y = 1- window(@gausswin,dimy,GW);
[wx wy] = meshgrid(w1y,w1x);
wg = fftshift(wx.*wy);
wg(1,1) = 1;
% Butterworth filter
Dim_XY = size(im);
CutoffFrequencies_XY = Dim_XY*0.99;
FilterOrder = 200;
wb = FilterButterworth(Dim_XY,CutoffFrequencies_XY,FilterOrder);
wb(1,1) = 1;
w = wg;
%% Read flat fields
for nn = 3:-1:1
    flat = FilterHotPixel(double(imread(sprintf('%sproj_%05u.tif',DataPath,1360+nn+1))),HPT);
    flatfilt(:,:,nn) = real(ifft2(fft2(flat).*w));
end
%% Fourier transform
imft   = fft2(im);
imfilt  = real(ifft2(imft.*w));
%domain(im);domain(imfilt);domain(flat(:,:,1));domain(flatfilt(:,:,1))
int1 = imfilt./flatfilt(:,:,1);
%ishow(imfilt);ishow(flatfilt(:,:,1));ishow(int1)
fprintf('Computing median: ');tic
flat = median(flatfilt,3);
fprintf('%f s.\n',toc)
int = imfilt./flat;
domain(im);domain(imft);domain(int);domain(flat)
ishow(flat);ishow(int)