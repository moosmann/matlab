function  [phi0,phi1,phi2] = rec(dataslice,padding,alpha,padvalue,renormalize);
% Phase retrieval algorithm to first and second order in distance z. Input
% is the intensity contrast at given distance. Output is the phase map of
% the first order (Bronnikov) and two correction terms to the first order
% result. Input data is zero-padded to the next power of 2 (or bigger if
% desired), and clipped to original size after finishing the phase
% retrieval. The algorithm needs a regularization due to the singularity
% encountered when inverting the Laplacian. This effects the absolut value
% of the retrieved phase. Therefore the result is renormalized to the
% intervall [-1,0], and the corrections acoordingly adjusted.

if (nargin<4) || isempty(padvalue)
    padvalue = 0;
end;
if (nargin<5) || isempty(renormalize)
    renormalize = 1;
end;

% Dimensions of input data array.
[dim2,dim1] = size(dataslice);
dimx        = padding*2^nextpow2(dim1);
dimy        = padding*2^nextpow2(dim2);
% Filters.
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi         = fftshift(xi);
eta        = fftshift(eta);
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + 10^-alpha);
% Define intensity contrast function g=I(x,y,z)/I(x,y,0)-1.
g          = (dataslice-mean(mean(dataslice)));
% Zero-pad intensity contrast
g          = padarray(g,[(dimy-dim2)/2,(dimx-dim1)/2],padvalue,'both');
% Fourier transform of g
g_ft       = fft2(g,dimy,dimx);
% Compute FT[phi] by inverting the laplacian.
phi0_ft    = 1/(2*pi)*inv_lap.*g_ft;
% First order.
phi0       = ifft2(phi0_ft);
% Second order corrections.
iftxiphi_ft  = ifft2( xi.*phi0_ft,dimy,dimx);
iftetaphi_ft = ifft2(eta.*phi0_ft,dimy,dimx);
phi1         = ifft2(inv_lap.*fft2( ...
               -1/(4*pi)*g.^2 ...
               -1/2* iftxiphi_ft.*ifft2( xi.*g_ft) ...
               -1/2*iftetaphi_ft.*ifft2(eta.*g_ft) ...
                               ,dimy,dimx),dimy,dimx,'symmetric');
phi2         = -pi/2*(iftxiphi_ft.^2 + iftetaphi_ft.^2);
%% Iterating the result. Does not seem to work most probably due to the
%% compulsory renormalization.
%while iterations > 0
%phi           = normat(phi0 + phi1 + phi2);
%phi_ft        = fft2(phi);
%iftxiphi_ft   = ifft2( xi.*phi_ft,dimy,dimx);
%iftetaphi_ft  = ifft2(eta.*phi_ft,dimy,dimx);
%phi1          = -pi*ifft2(inv_lap.*fft2( ...
%                 ifft2(lap.*phi_ft).^2 ... 
%               + iftxiphi_ft.*ifft2( xi.*lap.*phi_ft) ...
%               +iftetaphi_ft.*ifft2(eta.*lap.*phi_ft)     ),dimy,dimx);
%phi2          = -pi/2*(iftxiphi_ft.^2 + iftetaphi_ft.^2);
%iterations    = iterations - 1;
%end;
% Take real part and clip zero-padded matrices to original size
if dimx~=dim1 | dimy~=dim2
fprintf(['Zero-padded %u x %u data to %u x %u and clipped to original ' ...
         'size.\n'],dim1,dim2,dimx,dimy);
end;
ycut = 1+dimy/2-dim2/2:dimy/2+dim2/2;
xcut = 1+dimx/2-dim1/2:dimx/2+dim1/2;
phi0      = real(phi0(ycut,xcut));
phi1      = real(phi1(ycut,xcut));
phi2      = real(phi2(ycut,xcut));
% Renormalize phase map to intervall [-1,0]
if renormalize
minval = min(min(phi0+phi1+phi2));
maxval = max(max(phi0+phi1+phi2));
phi0   = (phi0-maxval)/(maxval-minval);
phi1   = (phi1)/(maxval-minval);
phi2   = (phi2)/(maxval-minval);
end;
