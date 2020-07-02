function erplot(sli,alphamax,pts,cut1,cut2);
%print:sum_i |phi0_i|(alpha),sum_i |phi_i|(alpha),sum_i (|phi0_i|-|phi_i|)(alpha)

for ii=0:pts;
  x(ii+1)=ii/pts*alphamax;
  [med0(ii+1),med1(ii+1),diff0(ii+1),diff(ii+1)]=recbron(sli,1,x(ii+1));
end;

%figure: median0,medan1
int=1:pts+1;
figure,plot(x(int),med0(int),'blue',x(int),med1(int),'red');
%figure: median0-median1
int=1:pts+1;
figure,plot(x(int),med0(int)-med1(int),'blue');
%figure: zooming in of figure above
if cut1~=0 & cut2~=0
int=1+(floor(cut1*pts/alphamax):floor(cut2*pts/alphamax));
figure,plot(x(int),med0(int)-med1(int),'blue');
end;


%figure: diff(phi_exact-phi0),diff(phi_exact-phi)
int=1:pts+1;
figure,plot(x(int),diff0(int),'blue',x(int),diff(int),'red');
%figure: zooming in of figure above
if cut1~=0 & cut2~=0
int=1+(floor(cut1*pts/alphamax):floor(cut2*pts/alphamax));
figure,plot(x(int),er1(int),'blue',x(int),er2(int),'red');
end;



