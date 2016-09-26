function [I,fprop,fu,fph,xi,eta] = Propagate_Grating(distance);
% Compute intensity pattern at distance z in Fresnel theory for
% monochromatic incident plane wave of wave length lambda.

if (nargin<2) || isempty(lambda)
    %    lambda = 0.886e-10;
    lambda = 0.5166e-10;
end;
if (nargin<1) || isempty(distance)
    distance = 1.0;
end;

% PARAMETERS.
% Detector resolution.
dx   = 2e-6;
fprintf(1,'distance=%g, lambda=%g, lambda*z/dx^2=%g\n',distance,lambda, ...
        lambda*distance/dx^2);
% Define object function to propagate.
% grating 2.4e-6, sampling points 
px = 2.4e-6;
sx = 128;
object      = Grating([px,sx],[px,sx],0.496e-10,[3.915e-7,0]);
[dimx,dimy] = size(object);
padx        = 1*dimx;
pady        = padx;
object      = padarray(object,[(padx-dimx)/2,(pady-dimy)/2],0,'both');
% Propagation.
[xi,eta] = meshgrid(-1/2:1/padx:1/2-1/padx,-1/2:1/pady:1/2-1/pady);
xi       = fftshift(xi);
eta      = fftshift(eta);
fprop    = (exp(-i*pi*lambda*distance/(dx^2)*(xi.^2+eta.^2)));
fu       = fft2(object);
% Intensity pattern.
I     = abs((ifft2(fprop.*fu))).^2;
xcut  = 1+(padx-dimx)/2:(padx+dimx)/2;
ycut  = 1+(pady-dimy)/2:(pady+dimy)/2;
I     = I(xcut,ycut);
domain(I);
ishow(I);