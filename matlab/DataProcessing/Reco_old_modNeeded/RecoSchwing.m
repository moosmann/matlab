function [phi] ...
    = RecoSchwing(dataslice,sn_smax,lambda,distance,pixelsize,padding,padvalue,compute_correction);

% Phase retrieval algorithm to leading, next-to-leading and
% next-to-next-to-leading order in z (sample-detector distance). Input is
% pure phase-contrast intensity pattern given at z. Three output stacks:
% First stack includes the phase map of the leading order (Bronnikov) and
% next-to-leading. Second stack is the first stack before inverting the
% Laplacian. Input data is zero-padded to the next power of 2 (or bigger if
% desired), and clipped to original size after finishing the phase
% retrieval. The algorithm needs a regularization due to the singularity
% encountered when inverting the Laplacian. This effects the absolut value
% of the retrieved phase.

if nargin<3,lambda=1;end
if nargin<4,distance=1;end
if nargin<5,pixelsize=1;end
if nargin<6,padding = 1;end
if nargin<7,padvalue = 0;end
if nargin<8,compute_correction=1;end

prefactor = pixelsize^2/(2*pi*lambda*distance);
sn    = sn_smax(1);
smax  = sn_smax(2);
ds    = smax/sn;
% Dimensions of input data array.
[dim1,dim2] = size(dataslice);
dimx        = padding*2^nextpow2(dim1);
dimy        = padding*2^nextpow2(dim2);
%fprintf('initial resolution: %g x %g, padded resolution: %g x %g\n',dim1,dim2,dimx,dimy);
% Filters.
[xi,eta,s]   = meshgrid(-1/2:1/(dimy):1/2-1/(dimy),-1/2:1/(dimx):1/2-1/(dimx),smax*(1:sn)/sn);
% xi = matrix of identic row ranging accroding to the first entry of
% meshgrid, the rows are repeated according to the length of the vector
% of the second entry of meshgrid. eta analague to xi.
xi         = fftshift(xi);
eta        = fftshift(eta);
xieta      = xi.*eta;
lap        = xi.^2 + eta.^2;
selap        = s.*exp(-s.*sqrt(lap));
% Define intensity contrast function g=I(x,y,z)/I(x,y,0)-1.
% For pure phase objects, <g> states the conservation of flux <I>=1.
g          = dataslice-mean(dataslice(:));
% Zero-pad intensity contrast
g          = padarray(g,[(dimx-dim1)/2,(dimy-dim2)/2],padvalue,'both');
% Fourier transform of g
g_fts = repmat(fft2(g),[1,1,sn]);
% Compute FT[phi] by inverting the laplacian.
%phi1_ft    = inv_lap.*g_ft;
% Leading order: Bronnikov.
phi1 = ds*sum(ifft2(selap.*g_fts),3);

phi1        = prefactor*real(phi1);
phi1        = phi1-mean(phi1(:));

if compute_correction==1
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
    %lap_phi2    = lap_phi2  - mean(lap_phi2(:));
    phi2        = ifft2(inv_lap.*fft2(lap_phi2));
    
    phi2        = prefactor*real(phi2);
    phi2        = phi2-mean(phi2(:));
    % Stack Images.
    phi         = cat(3,phi1,phi2);
    
else
    phi = cat(3,phi1,phi1);
end

% Take real part and clip zero-padded matrices to original size
if dimx~=dim1 || dimy~=dim2
    xcut  = 1+ceil((dimx-dim1)/2):floor((dimx+dim1)/2);
    ycut  = 1+ceil((dimy-dim2)/2):floor((dimy+dim2)/2);
    %fprintf('xcut min: %g, xcut max: %g, range: %g \n',min(xcut),max(xcut),max(xcut)-min(xcut));
    %fprintf('ycut min: %g, ycut max: %g, range: %g \n',min(ycut),max(ycut),max(ycut)-min(ycut));
    phi  = real(phi(xcut,ycut,:));
    phi(:,:,1) = phi(:,:,1) - mean(mean(phi(:,:,1)));
    phi(:,:,2) = phi(:,:,2) - mean(mean(phi(:,:,2)));
else
    phi  = real(phi);
end


