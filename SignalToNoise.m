ca;
dimx = 1024;
dimy = dimx;
phamax = 0.3;
PreNoise = 00000;
[x y] = meshgrid(-1/2:1/dimx:1/2-1/dimx);
w1 = window(@gausswin,dimx,6);
w1 = 1 - fftshift(w1);
[wx wy] = meshgrid(w1);
w = wx.*wy;
domain(w)
%ishow([wx wy w])
sinx = sin(pi*2*x);
domain(sinx)
domain(sign(sinx))
pha1 = 1+phamax/2*sinx.*w;
pha2 = 1+phamax/2*imfilter(sign(sinx),fspecial('gaussian',[12 12],2)).*w;
domain(pha1),domain(pha2)
ishow([pha1 pha2])
edp = [20 0.5 2.2e-6];
if PreNoise(1)
    nmap = 1/PreNoise*double(imnoise(PreNoise*ones(dimx,'uint16'),'poisson'));
    int1 = Propagation2(pha1,nmap,edp,1,'symmetric');
    int2 = Propagation2(pha2,nmap,edp,1,'symmetric');
else
    int1 = Propagation(pha1,edp,2,'symmetric');
    int2 = Propagation(pha2,edp,2,'symmetric');
end
domain(int1);domain(int2)
ishow(int1),ishow(int2)

% im = ones(dimx);
% NoiseLevel = 1000;
% im2 = uint16(NoiseLevel*im);
% 
% impn = double(imnoise(im2,'poisson'));
% imgn = double(imnoise(im2,'gaussian',0,0.000000233));

%int = 1/NoiseLevel*double(imnoise(uint16(NoiseLevel*int),'poisson'));