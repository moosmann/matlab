function routine_phan(padding,alpha);

cd ~/data/phantom;
load head;

%dimensions, number of projections
dim1 = 256;
dim2 = 256;
nop  = 360;

dimx            = padding*dim1;
dimy            = padding*dim2;

%read phantom data
imstack    = double(mexVolRead('proj40cm',[dim1 dim2 nop],'float32'));

%filters:
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi         = fftshift(xi);
eta        = fftshift(eta);
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + 10^-alpha);
    
t1 = tic;
    
%looping the stack
for ii=1:nop
    
%zero-pad image and define g
g       = padarray(imstack(:,:,ii)-1,[(dimy-dim2)/2,(dimx-dim1)/2],0,'both');

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


cd /media/tomoraid/moosmann/phantom_bronnikov;
pmedf_write(sprintf('phantom-phase-bronnikov_%4u.edf',1000+ii),head,phi0);
pmedf_write(sprintf('phantom-phase-corrections_%4u.edf',1000+ii),head,phi1);
pmedf_write(sprintf('phantom-phase-bron-plus-corr_%4u.edf',1000+ii),head,phi0+phi1);


end;

t2 = toc(t1);
t2/nop