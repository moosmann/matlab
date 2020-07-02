function [phase,int,lo,nlo,ldp,phialpha]=LenaSchwingBRO(PhotonCountsForAddingNoise,compute_correction,blurring,sorder_smaxfac)
 
%[phase,int,lo,nlo,ldp]=lena(PhotonCountsForAddingNoise,compute_correction,blurring)
    
    if nargin<1,PhotonCountsForAddingNoise=0;end;
    if nargin<2,compute_correction = 0;end;
    if nargin<3,blurring = [0 0 0];end;
    if nargin<4,sorder_smaxfac = [150 32];end;
% Parameters.
distance    = .5;
energy      = 30; %keV
lambda      = EnergyConverter(energy);
pixelsize   = .6e-6; %m
paddingProp = 4;
filestring  = '/home/moosmann/lena.tif';
alpha       = 12;
paddingReco = 1;
padvalue    = 0; 
ldp         = [lambda distance pixelsize];
sorder      = sorder_smaxfac(1);
smaxfac     = sorder_smaxfac(2);
resolution  = 1024;

pic = 6*normat(double(imread(filestring)));
[dimx dimy] = size(pic);
phase = zeros(2*dimx,2*dimy);
phase(257:768,257:768) = pic; 
if blurring(1)>0,
    hsizex = blurring(1);
    hsizey = blurring(2);
    sigma  = blurring(3);
    phase  = imfilter(phase,fspecial('gaussian',[hsizex hsizey],sigma));
end;
%phase = imfilter(2*SiemensStar(2048,256),fspecial('gaussian',[8 8],1.5));
int   = Propagation(phase,distance,lambda,pixelsize,paddingProp,0);
phase = phase - mean(phase(:));
if PhotonCountsForAddingNoise>0,
    counts = PhotonCountsForAddingNoise;
    int = 1/counts*double(imnoise(uint16(counts*int),'poisson'));
    Domain(int);
end;

%phi =%Reco(int,alpha,lambda,distance,pixelsize,paddingReco,padvalue,compute_correction);
phi = RecoSchwingGauss(int,[sorder floor(smaxfac*resolution)],lambda,distance,pixelsize,paddingReco,padvalue,compute_correction);

lo  = phi(:,:,1);
nlo = phi(:,:,2);
phialpha = Reco(int,alpha,lambda,distance,pixelsize,1,0,0,1);
