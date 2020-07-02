function schwing(dat,sli,pad,sn,smax)
%phase retrieval to first and second order

ds    = smax/sn;

%dimensions of input data array
[dimy,dimx,~] = size(dat);
dimx  = pad*dimx;
dimy  = pad*dimy;

%pick 2D slice of 3D array stack of projections
g          = dat(:,:,sli)-1;

%filters:
[xi,eta,s]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy),smax*(1:sn)/sn);
xi           = fftshift(fftshift( xi,1),2);
eta          = fftshift(fftshift(eta,1),2);
selap        = s.*exp(-s.*sqrt(xi.^2+eta.^2));

%FT of g
g_ft  = fft2(g,dimy,dimx);
g_fts = repmat(g_ft,[1,1,sn]);

%zeroth order result
phi0  = 1/(2*pi)*ds*sum((ifft2(selap.*g_fts)),3);

phi0  = real(phi0);

%first order result: corrections
phi11 = -1/(4*pi)*ds*sum((ifft2(selap.*repmat(fft2(g.^2,dimy,dimx),[1,1,sn]))),3);
phi12 = -1/(4*pi)*ds^3*sum(real(ifft2(selap.*repmat(fft2( ...
    sum((ifft2( xi.*selap.*g_fts)),3).*ifft2( xi(:,:,1).*g_ft) ...
    +sum((ifft2(eta.*selap.*g_fts)),3).*ifft2(eta(:,:,1).*g_ft) ...
    ),[1,1,sn]))),3);
phi13 = -1/(8*pi)*ds^2*(sum((ifft2(xi.*selap.*g_fts)),3).^2 +sum((ifft2(eta.*selap.*g_fts)),3).^2);

phi11 = real(phi11);
phi13 = real(phi13);

%clipping
if pad>1
    phi0  = phi0(1:dimy/pad,1:dimx/pad);
    phi11 = phi11(1:dimy/pad,1:dimx/pad);
    phi12 = phi12(1:dimy/pad,1:dimx/pad);
    phi13 = phi13(1:dimy/pad,1:dimx/pad);
end

%figures
if 0
    figure,imshow([phi0 ],[],'InitialMagnification','fit'),colorbar;
    figure,imshow([phi11],[],'InitialMagnification','fit'),colorbar;
    figure,imshow([phi12],[],'InitialMagnification','fit'),colorbar;
    figure,imshow([phi11+phi12],[],'InitialMagnification','fit'),colorbar;
    figure,imshow([phi13],[],'InitialMagnification','fit'),colorbar;
end


