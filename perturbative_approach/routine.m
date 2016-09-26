function routine(first,last,padding,alpha);

cd /media/tomoraid/rolo/IHR_2009_06/hs_tomo/binsengras_2_pt1/binsengras_2_marker_p1/optfilt;
name = 'optfilt_';

%get dimensions
[im.header,im.data] = pmedf_read(sprintf('%s%4u.edf',name,first));
[dim2,dim1] = size(im.data);
dimx        = padding*dim1;
dimy        = padding*dim2;

%filters:
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi         = fftshift(xi);
eta        = fftshift(eta);
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + 10^-alpha);
    
t1 = tic;
    
for ii=first:last
    
%get filename string
im_name = sprintf('%s%4u.edf',name,ii);

%read data and header
[im.header,im.data] = pmedf_read(im_name);
im.data = im.data - mean(mean(im.data));
%zero-pad image and define g
n       = mean(mean(im.data(950:980,:))); 
im.data(:,930:end) = n;
g       = padarray(im.data,[(dimy-dim2)/2,(dimx-dim1)/2],0,'both');

%FT of g
g_ft    = fft2(g,dimy,dimx);

%compute FT[phi] by inverting the  laplacian
phi0_ft = 1/(2*pi)*inv_lap.*g_ft;

%first order
phi0    = ifft2(phi0_ft,dimy,dimx);

%second order corrections:
iftxiphi_ft  = ifft2( xi.*phi0_ft,dimy,dimx);
iftetaphi_ft = ifft2(eta.*phi0_ft,dimy,dimx);

phi1         =  ifft2(inv_lap.*fft2( ...
               -1/(4*pi)*g.^2 ... 
               -1/2* iftxiphi_ft.*ifft2( xi.*g_ft) ...
               -1/2*iftetaphi_ft.*ifft2(eta.*g_ft) ...
                               ,dimy,dimx),dimy,dimx,'symmetric') ...
               -pi/2*(iftxiphi_ft.^2 + iftetaphi_ft.^2);

%take real part and clip zero-padded matrices to original size
if padding > 1
phi0   = real(phi0(1+dim2/2:3*dim2/2,1+dim1/2:3*dim1/2));
phi1   = real(phi1(1+dim2/2:3*dim2/2,1+dim1/2:3*dim1/2));
else
phi0   = real(phi0);
phi1   = real(phi1);
end;    

%renormalize phase
phimin = min(min(phi0+phi1));
phimax = max(max(phi0+phi1));
phi0   = (phi0-phimin)/(phimax-phimin);
phi1   = (phi1)/(phimax-phimin);

cd ../optfilt_phase_reco_norm2
pmedf_write(sprintf('bro_%4u.edf',ii),im.header,phi0);
pmedf_write(sprintf('brocor_%4u.edf',ii),im.header,phi0+phi1);
pmedf_write(sprintf('cor_%4u.edf',ii),im.header,phi1);
cd ../optfilt;


end;

t2 = toc(t1);
fprintf(1,'Total time elapsed: %u min\n',t2/60);
