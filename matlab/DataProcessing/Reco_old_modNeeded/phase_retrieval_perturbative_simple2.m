function  [phi0,phi1,phi2,g] = phase_retrieval_perturbative_simple2(dataslice,padding,alpha,beta);
%phase retrieval algorithm to first and second order in distance z


%dimensions of input data array
[dim2,dim1] = size(dataslice);
dimx        = padding*2^nextpow2(dim1);
dimy        = padding*2^nextpow2(dim2);

%filters:
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi         = fftshift(xi);
eta        = fftshift(eta);
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + 10^-alpha)+beta;
%detector psf
psf        = exp(-sqrt(2)*((xi/.5).^2+(eta/.5).^2));
ranges(psf);

%mean of data
datamean   = mean(mean(dataslice));

%zero-pad data
datapad    = padarray(dataslice,[(dimy-dim2)/2,(dimx-dim1)/2],0,'both');

%FT(g)/(point spread function)-FT(mean of data)
g_ft       = fft2(datapad,dimy,dimx)./psf - fft2(datamean,dimy,dimx) ;

g          = ifft2(g_ft);

%compute FT[phi] by inverting the  laplacian
phi0_ft    = 1/(2*pi)*inv_lap.*g_ft;

%first order
phi0       = ifft2(phi0_ft);

%second order corrections:
iftxiphi_ft  = ifft2( xi.*phi0_ft,dimy,dimx);
iftetaphi_ft = ifft2(eta.*phi0_ft,dimy,dimx);

phi1        = ifft2(inv_lap.*fft2( ...
               -1/(4*pi)*g.^2 ...
               -1/2* iftxiphi_ft.*ifft2( xi.*g_ft) ...
               -1/2*iftetaphi_ft.*ifft2(eta.*g_ft) ...
                               ,dimy,dimx),dimy,dimx,'symmetric');
phi2         = -pi/2*(iftxiphi_ft.^2 + iftetaphi_ft.^2);

%%iterations
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

%take real part and clip zero-padded matrices to original size
ycut = 1+dimy/2-dim2/2:dimy/2+dim2/2;
xcut = 1+dimx/2-dim1/2:dimx/2+dim1/2;
phi0      = real(phi0(ycut,xcut));
phi1      = real(phi1(ycut,xcut));
phi2      = real(phi2(ycut,xcut));

%renormalize matrix
minval = min(min(phi0+phi1+phi2));
maxval = max(max(phi0+phi1+phi2));
phi0   = (phi0-maxval)/(maxval-minval);
phi1   = (phi1)/(maxval-minval);
phi2   = (phi2)/(maxval-minval);
