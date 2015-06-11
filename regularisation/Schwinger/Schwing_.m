function [lo,nlo,star,int,ldp,phialpha,xi,eta,s] = Schwing(sorder,smaxfac,compute_correction)

% Test script for Schwinger regularization procedure.
    
% Parameters.
energy = 30;%keV
lambda = EnergyConverter(energy); %m
distance = 0.3; %m
pixelsize = 0.7e-6; %m
resolution = 1024;
spokes = 32;
hsize = 8;
sigma = 1.5;
delta = 0.2e-7;
thickness = 0.001;
ldp = [lambda distance pixelsize];
alpha = 12;

% Phase object (Siemens star) construction.
star = SiemensStar(resolution,spokes);
star = imfilter(star,fspecial('gaussian',[hsize hsize],sigma));
star = 2*pi/lambda*delta*thickness*star;
% Propagation.
int = Propagation(star,distance,lambda,pixelsize,2,0);
star = star - mean(star(:));
% Phase retrieval.
tic;
[phi,xi,eta,s]=RecoSchwingGauss(int,[sorder floor(smaxfac*resolution)],lambda,distance,pixelsize,1,0,compute_correction);
treco = toc;
phialpha = Reco(int,alpha,lambda,distance,pixelsize,1,0,0,1);
fprintf(1,'Time for phase retrieval: %gs\n',treco);
lo  = phi(:,:,1);
nlo = phi(:,:,2);
er  = abs(star - lo);
mer = mean(er(:));
% Print domains.
Domain(star);
Domain(lo);
% Print figures.
%ishow(lo);
%ishow(er);
% Save figures.

fprintf(1,'Mean error per pixel: %.12g\n\n',mer);
