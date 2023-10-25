function [phi0,phi11,phi12,phi13] = phase_retrieval_perturbative_simple(data,slice,~,alpha,~)
%phase retrieval algorithm to first and second order in distance z


%dimensions of input data array
[dim2,dim1,~] = size(data);
dimx             = dim1;
dimy             = dim2;

%pick 2D slice of 3D array stack of projections
g          = data(:,:,slice)-1;

%filters:
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi         = fftshift(xi);
eta        = fftshift(eta);
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + 10^-alpha);

%zero-pad intensity contrast
g          = padarray(g,[(dimy-dim2)/2,(dimx-dim1)/2]);

%FT of g
g_ft       = fft2(g);

%compute FT[phi] by inverting the  laplacian
phi0_ft    = 1/(2*pi)*inv_lap.*g_ft;

%first order
phi0       = ifft2(phi0_ft);

%second order corrections:
phi11      = -1/(4*pi)*ifft2(inv_lap.*fft2(g.^2,dimy,dimx));
phi12      = -1/2*(ifft2(inv_lap.*fft2(ifft2( xi.*phi0_ft).*ifft2( xi.*g_ft) ...
				      +ifft2(eta.*phi0_ft).*ifft2(eta.*g_ft) ...
				       )));
phi13      = -pi/2*((ifft2( xi.*phi0_ft).^2 ... 
		   + ifft2(eta.*phi0_ft).^2));

phi        = phi0 + phi11 + phi12 + phi13;

%take real part and clip zero-padded matrices to original size
phi0       = real( phi0(1:dim2,1:dim1));
phi11      = real(phi11(1:dim2,1:dim1));
phi12      = real(phi12(1:dim2,1:dim1));
phi13      = real(phi13(1:dim2,1:dim1));
