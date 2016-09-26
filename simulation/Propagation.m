function [intensity,intensity_padded,fprop,fu,xi,eta] = ...
    Propagation(phase_object,EnergyDistancePixelsize,padding,MethodOrPadValue,PrintDomains) 

% Compute intensity pattern at distance z in Fresnel theory for
% monochromatic incident plane wave of given energy and assumming given pixelsize.
%
% [intensity,intensity_padded,fprop,fu,xi,eta] = 
% Propagation(phase_object,EnergyDistancePixelsize,padding,MethodOrPadValue,PrintDomains); 

if nargin<2
    EnergyDistancePixelsize(1)=25; % in keV
    EnergyDistancePixelsize(2)=1; % in m
    EnergyDistancePixelsize(3)= 0.36e-6;
end; % in m
if nargin<3
    padding = 1;
end;
if nargin<4
    MethodOrPadValue = 'symmetric';
end;
if nargin<5
    PrintDomains = 1;
end;

% Parameters.
energy      = EnergyDistancePixelsize(1); % in keV
lambda      = EnergyConverter(energy); % in m
distance    = EnergyDistancePixelsize(2); % in m
pixelsize   = EnergyDistancePixelsize(3); % in m
[dimx,dimy] = size(phase_object);

% Program.
% Pad object.
padx        = padding*dimx;
pady        = padding*dimy;
phase_object= padarray(phase_object,[(padx-dimx)/2,(pady-dimy)/2],MethodOrPadValue,'both');
% Fourier cooridnates.
[xi,eta] = meshgrid(-1/2:1/pady:1/2-1/pady,-1/2:1/padx:1/2-1/padx);
xi       = fftshift(xi);
eta      = fftshift(eta);
% Propagator.
fprop    = (exp(-1i*pi*lambda*distance/(pixelsize^2)*(xi.^2+eta.^2)));
fu       = fft2(exp(1i*phase_object));
% Intensity pattern.
intensity_padded     = abs((ifft2(fprop.*fu))).^2;
xcut  = 1+(padx-dimx)/2:(padx+dimx)/2;
ycut  = 1+(pady-dimy)/2:(pady+dimy)/2;
intensity     = intensity_padded(xcut,ycut);
% Print parameters and domains.
if PrintDomains,
    fprintf(1,'distance=%gm, lambda=%gm, pixelsize=%g, resolution=[%u,%u], lambda*distance/pixelsize^2=%g\n', ...
        distance,lambda,pixelsize,dimx,dimy,lambda*distance/pixelsize^2);
    domain(intensity,1,'intensity');
    %domain(phase_object,'exact phase');
    %domain(phase_object/2/phi,'phase/2pi');
end;


