function [phi0,phi1,x] = alphastacks(dataslice,padding,alphamax,pts);

for ii=0:pts;
  x(ii+1)                  = ii/pts*alphamax;
  [phi0(:,:,ii+1),b,c] = rec(dataslice,padding,x(ii+1),0);
  phi1(:,:,ii+1)=b+c;

end;
