function [f,s,xi,eta] = SineFilter(resolution,EnergyDistancePixelsize,width, shift)

% function [f,s,xi,eta] = SineFilter(resolution,EnergyDistancePixelsize,width, shift)    
    
    if nargin<1, resolution = [1024,1024];end;
    if nargin<2, EnergyDistancePixelsize(1)=30;EnergyDistancePixelsize(2)=1;EnergyDistancePixelsize(3)=1.1e-6;end;
    if nargin<3, width=0.2;end;
    if nargin < 4
        shift = 0;
    end

dimx      = resolution(1);
dimy      = resolution(2);
lambda    = EnergyConverter(EnergyDistancePixelsize(1));
distance  = EnergyDistancePixelsize(2);
pixelsize = EnergyDistancePixelsize(3);
[xi,eta]  = meshgrid(-1/2:1/dimy:1/2-1/dimy,-1/2:1/dimx:1/2-1/dimx);
args      = pi*lambda*distance/pixelsize^2*(xi.^2+eta.^2);
s         = sin(args - shift);
f         = ones(dimx,dimy);
f((s.^2<width)&args>pi/2) = 0;

