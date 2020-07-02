function  [mb,mbc,erb,erbc,a,b] = ervsalphabeta(stack,slice,padding,amax,apts,bmax,bpts);
%error measure vs regularization parameters alpha and beta

    global phase;
    pha = phase(:,:,slice);
    im  = stack(:,:,slice);
for ii=0:apts;
  a(ii+1)          = ii/apts*amax;
  for jj=0:bpts;
  b(jj+1)          = jj/bpts*bmax;
  [phi0,phi1,phi2] = rec(im,padding,a(ii+1),b(jj+1));
  mb(ii+1,jj+1)    = mean(mean(abs(normat(phi0))));
  mbc(ii+1,jj+1)   = mean(mean(abs(normat(phi0+phi1+phi2))));
  erb(ii+1,jj+1)   = mean(mean(abs(normat(phi0)-normat(pha))));
  erbc(ii+1,jj+1)  = mean(mean(abs(normat(phi0+phi1+phi2)-normat(pha))));
  end;
end;

