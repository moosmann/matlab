function schwing2(dat,sli,pad,sn,smax);
%phase retrieval to first and second order

%dimensions of input data array
dims  = size(dat);
dimx  = dims(2);
dimy  = dims(1);

%pick 2D slice of 3D array stack of projections
g          = dat(:,:,sli)-1;

%filters:
[xi,eta,s]   = meshgrid(-1/2:1/(pad*dimx):1/2-1/(pad*dimx),-1/2:1/(pad*dimy):1/2-1/(pad*dimy),smax*(1:sn)/sn);
xi           = fftshift(fftshift( xi,1),2);
eta          = fftshift(fftshift(eta,1),2);
selap        = s.*exp(-s.*sqrt(xi.^2+eta.^2));

%FT of g

g_fts        = repmat(fft2(g,pad*dimy,pad*dimx),[1,1,sn]);

%zeroth order result
%sphi0        = ifft2(selap.*g_ft);
phi0         = 1/(2*pi)*sum((ifft2(selap.*g_fts)),3)/sn;

phi0 = real(phi0);

%first order result
%first correction of three
phi11        = -1/(4*pi)*sum((ifft2(selap.*repmat(fft2(g.^2,pad*dimy,pad*dimx),[1,1,sn]))),3)/sn;
%second correction
phi12        = -1/(4*pi)*sum(real(ifft2(selap.*repmat(fft2( ...
                sum((ifft2( xi.*selap.*g_fts)),3)/sn.*ifft2( xi(:,:,1).*g_fts(:,:,1)) ...
               +sum((ifft2(eta.*selap.*g_fts)),3)/sn.*ifft2(eta(:,:,1).*g_fts(:,:,1)) ...
                                                           ),[1,1,sn]))),3)/sn;
phi13        = -1/(8*pi)*((sum((ifft2(xi.*selap.*g_fts)),3)/sn).^2 ... 
                         +(sum((ifft2(eta.*selap.*g_fts)),3)/sn).^2);
%figran(phi11);
%figran(phi12);
%figran(phi13);
phi11 = real(phi11);
phi13 = real(phi13);

%clip region of interest
if pad>1
phi0  = phi0(1:dimy,1:dimx);
phi11 = phi11(1:dimy,1:dimx);
phi12 = phi12(1:dimy,1:dimx);
phi13 = phi13(1:dimy,1:dimx);
end;

%figures
if 1
figure,imshow([phi0 ],[],'InitialMagnification','fit'),colorbar;
figure,imshow([phi11],[],'InitialMagnification','fit'),colorbar;
figure,imshow([phi12],[],'InitialMagnification','fit'),colorbar;
%figure,imshow([phi11+phi12],[],'InitialMagnification','fit'),colorbar;
figure,imshow([phi13],[],'InitialMagnification','fit'),colorbar;
end;


