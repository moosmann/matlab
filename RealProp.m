function [fu]=RealProp(lambda,distance);

if (nargin<1) || isempty(lambda)
    lambda = 0.1e-9;
end;
if (nargin<2) || isempty(distance)
    distance = 50e-2;
end;

% Read Bronnikov phantom.
cd ~/data/phantom
object = (double(mexVolRead('phase',[256 256 360],'float32')));
object = 1e10*object(:,:,1);
[dim,dim]=size(object);
domain(object);

% Laboratory parameters.
alpha = sqrt(2*pi*distance/lambda);
fprintf(1,'distance=%g, lambda=%g, alpha=%g\n',distance,lambda,alpha);
alpha = 1;

% Propagation.
[xi,eta]=meshgrid(-1/2:1/dim:1/2-1/dim,-1/2:1/dim:1/2-1/dim);
xi=xi;eta=eta;
fprop=exp(-i*pi*lambda*distance*(xi.^2+eta.^2));
I=abs(ifft2(fprop.*fft2(exp(i*object))))^2;
domain(I);
ishow(I);