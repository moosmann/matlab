function [I,fprop,fu,fph,xi,eta] = Propagation_Volume(dmax,Nz,lambda,pixelsize);
% Compute intensity pattern at distance z in Fresnel theory for
% monochromatic incident plane wave of wave length lambda.

if (nargin<3) || isempty(lambda)
    % Wave length.
    lambda = 1e-10; %12.4 keV
end;
if (nargin<4) || isempty(pixelsize)
    % Detector resolution.
    pixelsize = 1e-6;
end;


%fprintf(1,'distance=%g, lambda=%g, lambda*z/dx^2=%g\n',distance,lambda, ...
%       lambda*distance/dx^2);
% Object function to propagate.
% Read Bronnikov phantom.
cd ~/data/phantom;
object      = (double(mexVolRead('phase',[256 256 360],'float32')));
object      = -1e2*object(:,:,1);
object      = 1e8*object;%-min(object(:));
[dimx,dimy] = size(object);
domain(object);
% Padding.
padx        = 1*dimx;
pady        = padx;
xcut        = 1+(padx-dimx)/2:(padx+dimx)/2;
ycut        = 1+(pady-dimy)/2:(pady+dimy)/2;
object      = padarray(object,[(padx-dimx)/2,(pady-dimy)/2],0,'both');
% Fourier transform of transmission function.
fu          = fft2(exp(i*object));
% Fresnel propagation.
[xi,eta] = meshgrid(-1/2:1/padx:1/2-1/padx,-1/2:1/pady:1/2-1/pady);
xi       = fftshift(xi);
eta      = fftshift(eta);
I        = zeros(padx,pady,Nz);
for nn = (0:Nz),
    distance = nn^3/Nz^3*dmax;
    fprop    = exp(-i*pi*lambda*distance/(pixelsize^2)*(xi.^2+eta.^2));
    % Intensity pattern.
    I(:,:,nn+1) = abs(ifft2(fprop.*fu)).^2;
    I(:,:,nn+1) = I(xcut,ycut,nn+1);
    clear fprop;
end;
ishow([I(:,:,1),I(:,:,Nz)]);