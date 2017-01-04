function [phi2] ... 
= Correction(phase,alpha,lambda,distance,pixelsize);
                                                  
if (nargin<3),alpha=12;end;
if (nargin<3),lambda=1;end;
if (nargin<4),distance=1;end;
if (nargin<5),pixelsize=1;end;

padding = 1;
padvalue = 0;
prefactor = pixelsize^2/(2*pi*lambda*distance);

% Dimensions of input data array.
[dim1,dim2] = size(phase);
dimx        = padding*2^nextpow2(dim1);
dimy        = padding*2^nextpow2(dim2);
% Filters.
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi         = fftshift(xi);
eta        = fftshift(eta);
xieta      = xi.*eta;
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + 10^-alpha);
% Compute correction.
phase      = phase/prefactor;
phi1       = phase-mean(phase(:));
phi1_ft    = fft2(phase);
% Next-to-leading order.
phi1_ftx    = fft(phi1,[],1);
phi1_fty    = fft(phi1,[],2);
phi_dx1     = ifft(eta   .*phi1_ftx,[],1);
phi_dx2     = ifft(eta.^2.*phi1_ftx,[],1);
phi_dx3     = ifft(eta.^3.*phi1_ftx,[],1);
phi_dy1     = ifft( xi   .*phi1_fty,[],2);
phi_dy2     = ifft( xi.^2.*phi1_fty,[],2);
phi_dy3     = ifft( xi.^3.*phi1_fty,[],2);
phi_dx1dy1  = ifft2(xieta.*phi1_ft);
phi_dx1dy2  = ifft2( xi.*xieta.*phi1_ft);
phi_dx2dy1  = ifft2(eta.*xieta.*phi1_ft);
lap_phi2    = -(phi_dx1.*(phi_dx3+phi_dx1dy2) + phi_dy1.*(phi_dx2dy1+phi_dy3) ...
            + phi_dx2.^2 + phi_dx2.*phi_dy2 + phi_dx1dy1.^2 + phi_dy2.^2);
lap_phi2    = real(lap_phi2);
phi2        = ifft2(inv_lap.*fft2(lap_phi2));
phi2        = prefactor*(phi2 - mean(phi2(:)));


% Take real part and clip zero-padded matrices to original size
if dimx~=dim1 || dimy~=dim2
    xcut  = 1+dimx/2-dim1/2:dimx/2+dim1/2;
    ycut  = 1+dimy/2-dim2/2:dimy/2+dim2/2;
    phi2  = real(phi2(xcut,ycut,:));
else
    phi2  = real(phi2);
end;
