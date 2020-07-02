function [rm,x,y]=ConcentricRings(dim,nyquistx)
% Create concentric rings of decreasing amplitude and period.
 
if (nargin<1) || isempty(dim)
    dim = 1024;
end;
if (nargin<2) || isempty(nyquistx)
    nyquistx = floor(dim/2);
end;

nyquisty = nyquistx;
dimx=dim;
dimy=dim;
dimhalf = ceil(dim/2);
% Radius matrix.
[x,y] = meshgrid(-dimhalf:dimhalf-1,-dimhalf:dimhalf-1);
r     = sqrt(x.^2+y.^2);
% 1D-stripes in fourier space.
if 0
st    = zeros(dimx,dimy);
st(dimhalf-nyquistx:dimhalf+nyquistx,:)=1;
ishow(st);
iftst=abs((fft(st)));
domain(iftst);
figure,plot(iftst(1:dimhalf,dimhalf))
end;
rm=1/dimx*sin(pi*nyquistx/dimx*x)./sin(pi/dimx*x);
rm(:,dimhalf+1)=nyquistx/dimx;
domain(rm);
l=rm(dimhalf,:);
figure,plot(l);
lf = fft(l);
domain(lf);
figure,plot(fftshift(real(lf)));


