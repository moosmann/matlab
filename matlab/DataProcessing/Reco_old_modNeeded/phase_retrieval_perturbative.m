function [phi,fts,lphi] = phase_retrieval_perturbative(dataslice,alpha,padding,padvalue,renormalize)
% Phase retrieval algorithm to leading, next-to-leading and
% next-to-next-to-leading order in z (sample-detector distance). Input is
% pure phase-contrast intensity pattern given at z. Output stack includes
% the phase map of the leading order (Bronnikov) (1.stack) and correction
% terms for next-to-leading (2.&3.stack) and next-to-next-to leading
% order(4.-12.stack). Input data is zero-padded to the next power of 2 (or
% bigger if desired), and clipped to original size after finishing the phase
% retrieval. The algorithm needs a regularization due to the singularity
% encountered when inverting the Laplacian. This effects the absolut value
% of the retrieved phase. Therefore the result is renormalized to the
% intervall [-1,0], and the corrections acoordingly adjusted.

if (nargin<4)||isempty(padvalue),padvalue = 0;end
if (nargin<5)||isempty(renormalize),renormalize = 0;end

% Dimensions of input data array.
[dim2,dim1] = size(dataslice);
dimx        = padding*2^nextpow2(dim1);
dimy        = padding*2^nextpow2(dim2);
phi         = zeros(dimx,dimy,12);
% Filters.
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi         = fftshift(xi);
eta        = fftshift(eta);
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + 10^-alpha);
% Define intensity contrast function g=I(x,y,z)/I(x,y,0)-1.
g          = dataslice-mean(mean(dataslice));
% Zero-pad intensity contrast
g          = padarray(g,[(dimy-dim2)/2,(dimx-dim1)/2],padvalue,'both');
% Fourier transform of g
g_ft       = fft2(g,dimy,dimx);
fts(:,:,1) = g_ft;
% Compute FT[phi] by inverting the laplacian.
phi0_ft    = inv_lap.*g_ft;
fts(:,:,2) = phi0_ft;
% Leading order: Bronnikov.
phi0         = ifft2(phi0_ft);
phi(:,:,1)   = phi0;
% Next-to-leading order.
iftxiphi_ft  = ifft2( xi.*phi0_ft,dimy,dimx);
iftetaphi_ft = ifft2(eta.*phi0_ft,dimy,dimx);
%iftxig_ft    = ifft2( xi.*g_ft);iftetag_ft   = ifft2(eta.*g_ft);
iftxig_ft    = ifft2( xi.*lap.*phi0_ft);iftetag_ft   = ifft2(eta.*lap.*phi0_ft);
lphi(:,:,1)  = -1/2* iftxiphi_ft.* iftxig_ft-1/2*iftetaphi_ft.*iftetag_ft;
phi(:,:,2)   = ifft2(inv_lap.*fft2(lphi(:,:,1)));
%lphi(:,:,2)  = -1/2*g.^2;
%lphi(:,:,2)  = -1/2*g.*ifft2(lap.*phi0_ft);
lphi(:,:,2)  = -1/2*ifft2(lap.*phi0_ft).^2;
phi(:,:,3)   = ifft2(inv_lap.*fft2(lphi(:,:,2)));
phi(:,:,4)   = -1/4*(iftxiphi_ft.^2 + iftetaphi_ft.^2);
phi22        = phi(:,:,3);
% Next-to-next-to-leading order.
phi(:,:,5)  = ifft2(-1/24*lap.*g_ft);
phi(:,:,6)  = -1/12*(iftxiphi_ft.*ifft2( xi.*fft2(-4*phi22)) ...
                 + iftetaphi_ft.*ifft2(eta.*fft2(-4*phi22)));
phi(:,:,7)  = ifft2(inv_lap.*fft2(- 1/6*g.^3));
phi(:,:,8)  = ifft2(inv_lap.*fft2(-1/12*g.*(iftxiphi_ft.*iftxig_ft ...
                                         +iftetaphi_ft.*iftetag_ft)));
phi(:,:,9)  = ifft2(inv_lap.*fft2(- 1/4*g.*(ifft2(lap.*fft2(-4*phi22)))));
phi(:,:,10) = ifft2(inv_lap.*fft2(- 1/3*(iftxig_ft.*ifft2( xi.*fft2(-4*phi22)) ...
                                      +iftetag_ft.*ifft2(eta.*fft2(-4*phi22)))));
phi(:,:,11) = ifft2(inv_lap.*fft2(- 1/6*(iftxiphi_ft.*ifft2( xi.*fft2(g.^2)) ...
                                       +iftetaphi_ft.*ifft2(eta.*fft2(g.^2)))));
phi(:,:,12) = ifft2(inv_lap.*fft2(- 1/6*(iftxiphi_ft.*ifft2( xi.*fft2( iftxig_ft.* iftxiphi_ft)) ...
                                       +iftetaphi_ft.*ifft2(eta.*fft2(iftetag_ft.*iftetaphi_ft)))));
phi(:,:,13) = ifft2(inv_lap.*fft2(-1/12*(iftxiphi_ft.*ifft2(lap.*xi.* fft2(-4*phi22)) ...
                                       +iftetaphi_ft.*ifft2(lap.*eta.*fft2(-4*phi22)))));
% Substract mean.
for ii=1:size(phi,3),phiii=phi(:,:,ii);phi(:,:,ii)=phiii-mean(phiii(:)); end
    clear phiii;

% Take real part and clip zero-padded matrices to original size
if dimx~=dim1 || dimy~=dim2
fprintf(['Zero-padded %u x %u data to %u x %u and clipped to original ' ...
         'size.\n'],dim1,dim2,dimx,dimy);
end
ycut  = 1+dimy/2-dim2/2:dimy/2+dim2/2;
xcut  = 1+dimx/2-dim1/2:dimx/2+dim1/2;
phi  = real(phi(ycut,xcut,:));
fts  = real(fts(ycut,xcut,:));
% Renormalize phase map to intervall [-1,0]
if renormalize
minval = min(phi(:));
maxval = max(phi(:));
phi(:,:,1) = phi(:,:,1)-maxval;
phi  = phi./(maxval-minval);
end
