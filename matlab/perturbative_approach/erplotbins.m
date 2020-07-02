function erplotbins(alphamax,pts,cut1,cut2);


for ii=0:pts;
  x(ii+1)=ii/pts*alphamax;
  [er1(ii+1),er2(ii+1)]=recbins(0,1,x(ii+1));
end;

int=1:pts+1;
figure,plot(x(int),er1(int),'blue',x(int),er2(int),'red');
if cut1~=0 & cut2~=0
int=1+(floor(cut1*pts/alphamax):floor(cut2*pts/alphamax));
figure,plot(x(int),er1(int),'blue',x(int),er2(int),'red');
end;
