function [phase,int,phiCTF,phiBRO,ldp,namestring,intFT]=LenaCTF(PhotonCountsForPoissonNoise,compute_correction,blurring,alpha,MaxPhaseShift,EnergyDistancePixelsize)

% [phase,int,phiCTF,phiBRO,ldp]=LenaCTF(PhotonCountsForPoissonNoise,compute_correction,blurring,alpha,MaxPhaseShift,EnergyDistancePixelsize)

    
    if nargin<1,PhotonCountsForPoissonNoise=0;end;
    if nargin<2,compute_correction = 0;end;
    if nargin<3,blurring = [0 0 0];end;
    if nargin<4,alpha = 12;end;
    if nargin<5,MaxPhaseShift = 3;end;
    if nargin<6,EnergyDistancePixelsize(1)=30;EnergyDistancePixelsize(2)= ...
            1;EnergyDistancePixelsize(3)= 1.1e-6;end;
% Parameters.
energy      = EnergyDistancePixelsize(1); %keV
lambda      = EnergyConverter(energy);
distance    = EnergyDistancePixelsize(2);
pixelsize   = EnergyDistancePixelsize(3); %m
paddingProp = 4;
filestring  = '/home/moosmann/lena.tif';
paddingReco = 1;
padvalue    = 0; 
ldp         = [lambda distance pixelsize];
counts_min  = 0;
counts_max  = 0;
        

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
if PhotonCountsForPoissonNoise>0,
    counts = PhotonCountsForPoissonNoise;
    if PhotonCountsForPoissonNoise*max(int(:))<2^16,
        int = imnoise(uint16(counts*int),'poisson');
        counts_min = min(int(:));
        counts_max = max(int(:));
        int = 1/counts*double(int);
    else
        fprintf(1,'Warning: Counts exceed 2^16\n');
    end;
end;
intmin = min(int(:));
intmax = max(int(:));
fprintf(1,['phi_max: %4.2f, Energy: %2ukeV, Distance: %4.2fm, ' ...
           'Counts: %5u [%5u,%5u], Intensity: [%7.5f,%7.5f] (%.5f)\n'], ... 
        MaxPhaseShift,energy,distance,PhotonCountsForPoissonNoise, ... 
        counts_min,counts_max,intmin,intmax,intmax-intmin);

phiCTF = RecoCTF(int,alpha,lambda,distance,pixelsize,paddingReco,padvalue,compute_correction);
phiBRO = Reco(int,alpha,lambda,distance,pixelsize,1,0,0,1);
intFT  = laf(int);

% Print/save images: Exact phase. Intensity. Retrieved phase.
namestring  = sprintf('_%04.2f',MaxPhaseShift);
namestring  = regexprep(namestring,'\.','p');
%mexVolWrite(sprintf('%spha%s.vol'    ,folder,namestring),normat(phase),'float32');
namestring  = sprintf('%s_E%02u_z%04.2f_c%05u',namestring,energy,distance,PhotonCountsForPoissonNoise);
namestring  = regexprep(namestring,'\.','p');
%mexVolWrite(sprintf('%sphase_int_FTint_CTF_BRO%s.vol',folder,namestring),cat(3,phase,int,intFT,phiCTF,phiBRO),'float32');


