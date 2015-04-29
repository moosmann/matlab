function out = Filter3D(vol,alphaCTF_alphaTIE,EnergyDistancePixelsize)
% Phase retrieval in tomographically reconstructed volume

%% Default arguments
if nargin<2
    alphaCTF_alphaTIE = 2.5;
end
if nargin < 3
    EnergyDistancePixelsize = [30 0.620 2.2e-6];
end

%vol = double(vol);
%% Parameters
[dimx, dimy, dimz] = size(vol);
dimx = single(dimx);
dimy = single(dimy);
dimz = single(dimz);
% Energy, distance, pixel size, wave length
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
% Prefactor needed for TIE and CTF retrieval.
prefactor = Pixelsize^2/(2*pi*lambda*Distance);
% Regularization parameters.
if size(alphaCTF_alphaTIE,2) == 1
    alphaCTF = alphaCTF_alphaTIE;
    alphaTIE = alphaCTF - log10(prefactor);
else
    %alphaCTF = alphaCTF_alphaTIE(1);
    alphaTIE = alphaCTF_alphaTIE(2);
end
vol = vol - mean(vol(:));
% %% Fourier transformation
% tic;
% volf = fftn(vol);
% fprintf('Fourier transform of %u x %u x %u volume computed in %gs.\n',dimx,dimy,dimz,toc);
%% Meshgrid
tic;
%[x y z] = meshgrid(single(-1/2:1/dimy:1/2-1/dimy),single(-1/2:1/dimx:1/2-1/dimx),single(-1/2:1/dimz:1/2-1/dimz));
[x, y, z] = meshgrid(-1/2:1/dimy:1/2-1/dimy,-1/2:1/dimx:1/2-1/dimx,-1/2:1/dimz:1/2-1/dimz);
%[x y] = meshgrid(-1/2:1/dimy:1/2-1/dimy,-1/2:1/dimx:1/2-1/dimx);
x = fftshift(x);
y = fftshift(y);
z = fftshift(z);
fprintf('3D meshgrid computed in %gs.\n',toc);
%% Phase retrieval
% CenFilWid = 40;
% lff = (1-exp(-((dimx*x).^2+(dimy*y).^2+(dimz*z).^2)/(2*CenFilWid^2)));
%CroFilWid = 8;
tic;
%out = prefactor*real(ifft2(fftn(vol)./(x.^2+y.^2+z.^2+10^-alphaTIE).*(1-exp(-((dimx*x).^2+(dimy*y).^2+(dimz*z).^2)/(2*CenFilWid^2))).*(1-exp(-(dimx*x/CroFilWid).^2/2)).*(1-exp(-(dimy*y/CroFilWid).^2/2)).*(1-exp(-(dimz*z/CroFilWid).^2/2))));
%out = prefactor*real(ifftn(fftn(vol).*(1-exp(-((dimx*x).^2+(dimy*y).^2+(dimz*z).^2)/(2*CenFilWid^2))).*(1-exp(-(dimx*x/CroFilWid).^2/2)).*(1-exp(-(dimy*y/CroFilWid).^2/2)).*(1-exp(-(dimz*z/CroFilWid).^2/2))));
%out = prefactor*real(ifft2(fftn(vol)./(x.^2+y.^2+z.^2+10^-alphaTIE).*(1-exp(-((dimx*x).^2+(dimy*y).^2+(dimz*z).^2)/(2*CenFilWid^2)))));
out = prefactor*real(ifftn(fftn(vol)./(x.^2+y.^2+z.^2+10^-alphaTIE)));
%out = prefactor*real(ifft2(fftn(vol)./(x.^2+y.^2+z.^2+10^-alphaTIE).*(1-exp(-((dimx*x).^2+(dimy*y).^2+(dimz*z).^2)/(2*CenFilWid)))));
fprintf('Phase retrieval (FT,filter,iFT,real) done in %gs.\n',toc);
nimplay(out)