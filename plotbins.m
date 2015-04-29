function plotbins(pad,alphamax,pts,cut1,cut2);
%figures of BINSENGRAS: median0,median;median0_zoomed,median_zoomed;median0-median

for ii=0:pts;
  x(ii+1)=ii/pts*alphamax;
  [med0(ii+1),med(ii+1)]=recbins(0,pad,x(ii+1));
end;

%figure: median0,median
int=1:pts+1;
figure('Name','Binsengras: mean_bron(blue), mean_we(red)'),plot(x(int),med0(int),'blue',x(int),med(int),'red');
if cut1~=0 & cut2~=0
int=1+(floor(cut1*pts/alphamax):floor(cut2*pts/alphamax));
figure('Name','Binsengras: mean_bron(blue),mean_we(red),zoomed'),plot(x(int),med0(int),'blue',x(int),med(int),'red');
end;
%figure: median0-median
int=1:pts+1;
figure('Name','Binsengras: mean_bron - mean_we'),plot(x(int),med0(int)-med(int),'blue');
