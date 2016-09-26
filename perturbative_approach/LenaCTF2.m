function [phase,int,lo,nlo,ldp,phialpha]=LenaCTF2(PhotonCountsForAddingNoise,compute_correction,blurring,alpha,MaxPhaseShift)
 
%[phase,int,lo,nlo,ldp]=lena(PhotonCountsForAddingNoise,compute_correction,blurring)
    
    if nargin<1,PhotonCountsForAddingNoise=0;end;
    if nargin<2,compute_correction = 0;end;
    if nargin<3,blurring = [0 0 0];end;
    if nargin<4,alpha = 12;end;
    if nargin<5,MaxPhaseShift = 3;end;
% Parameters.
distance    = 1;
energy      = 30; %keV
lambda      = EnergyConverter(energy);
pixelsize   = .6e-6; %m
paddingProp = 4;
filestring  = '/home/moosmann/lena.tif';
paddingReco = 1;
padvalue    = 0; 
ldp         = [lambda distance pixelsize];
resolution  = 1024;

pic = MaxPhaseShift*normat(double(imread(filestring)));
[dimx dimy] = size(pic);
phase = zeros(2*dimx,2*dimy);
phase(257:768,257:768) = pic; 
if blurring(1)>0,
    hsizex = blurring(1);
    hsizey = blurring(2);
    sigma  = blurring(3);
    phase  = imfilter(phase,fspecial('gaussian',[hsizex hsizey],sigma));
end;
phase = MaxPhaseShift*normat(phase);
int   = Propagation(phase,distance,lambda,pixelsize,paddingProp,0);
phase = phase - mean(phase(:));
if PhotonCountsForAddingNoise>0,
    counts = PhotonCountsForAddingNoise;
    int = 1/counts*double(imnoise(uint16(counts*int),'poisson'));
    Domain(int);
end;

phi = RecoCTF2(int,alpha,lambda,distance,pixelsize,paddingReco,padvalue,compute_correction);
lo  = phi(:,:,1);
nlo = phi(:,:,2);
phialpha = Reco(int,alpha,lambda,distance,pixelsize,1,0,0,1);
