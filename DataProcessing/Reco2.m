function [phi,lap_phi,inv_lap] ... 
= Reco2(dataslice,alpha,lambda,distance,pixelsize,padding,padvalue,gwidth,k_crit);
                                                  
% Phase retrieval algorithm to leading, next-to-leading and
% next-to-next-to-leading order in z (sample-detector distance). Input is
% pure phase-contrast intensity pattern given at z. Three output stacks:
% First stack includes the phase map of the leading order (Bronnikov) and
% next-to-leading. Second stack is the first stack before inverting the
% Laplacian. Input data is zero-padded to the next power of 2 (or bigger if
% desired), and clipped to original size after finishing the phase
% retrieval. The algorithm needs a regularization due to the singularity
% encountered when inverting the Laplacian. This effects the absolut value
% of the retrieved phase. Therefore the result is renormalized to the
% intervall [-1,0], and the corrections acoordingly adjusted.

if (nargin<3),lambda=1;end;
if (nargin<4),distance=1;end;
if (nargin<5),pixelsize=1;end;
if (nargin<6)||isempty(padding),padding = 1;end;
if (nargin<7)||isempty(padvalue),padvalue = 0;end;

prefactor = pixelsize^2/(2*pi*lambda*distance);
% Dimensions of input data array.
[dim1,dim2] = size(dataslice);
dimx        = padding*2^nextpow2(dim1);
dimy        = padding*2^nextpow2(dim2);
% Filters.
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
% xi = matrix of identic row ranging accroding to the first entry of
% meshgrid, the rows are repeated according to the length of the vector
% of the second entry of meshgrid. eta analague to xi.
xi         = fftshift(xi);
eta        = fftshift(eta);
xieta      = xi.*eta;
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + 10^-alpha);
inv_lap    = inv_lap.*(1-1./(1+exp((sqrt(dimx*dimy*lap)-k_crit)/gwidth)));
% Define intensity contrast function g=I(x,y,z)/I(x,y,0)-1.
% For pure phase objects, <g> states the conservation of flux <I>=1.
g          = dataslice-mean(mean(dataslice));
% Zero-pad intensity contrast
g          = padarray(g,[(dimx-dim1)/2,(dimy-dim2)/2],padvalue,'both');
% Fourier transform of g
g_ft       = fft2(g,dimy,dimx);
% Compute FT[phi] by inverting the laplacian.
phi1_ft    = inv_lap.*g_ft;
% Leading order: Bronnikov.
phi1        = ifft2(phi1_ft);
phi1        = phi1-mean(phi1(:));
ishow(real(phi1)),colorbar;
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
lap_phi2    = lap_phi2 - mean(lap_phi2(:));
phi2        = ifft2(inv_lap.*fft2(lap_phi2));
% Image stacks.
phi         = cat(3,prefactor*(phi1-mean(phi1(:))),prefactor*(phi2-mean(phi2(:))));
lap_phi     = cat(3,g,lap_phi2);
% Take real part and clip zero-padded matrices to original size
if dimx~=dim1 | dimy~=dim2,
    xcut  = 1+dimx/2-dim1/2:dimx/2+dim1/2;
    ycut  = 1+dimy/2-dim2/2:dimy/2+dim2/2;
    phi  = real(phi(xcut,ycut,:));
    lphi = real(lap_phi(xcut,ycut,:));
else,
    phi  = real(phi);
    lphi = real(lap_phi);
end;
