function [phase,int,phiCTF,phiBRO,namestring,EnergyDistancePixelsize]=Xenopus(PhotonCountsForPoissonNoise,ComputeCorrection,Blurring,Alpha,EnergyDistancePixelsize)

% [phase,int,phiCTF,phiBRO,EnergyDistancePixelsize,namestring,intFT] =
% Xenopus(PhotonCountsForPoissonNoise,ComputeCorrection,Blurring,Alpha,EnergyDistancePixelsize) 

    
% DEFAULTS.
    if nargin<1,PhotonCountsForPoissonNoise=0;end;
    if nargin<2,ComputeCorrection = 0;end;
    if nargin<3,Blurring = [0 0 0];end;
    if nargin<4,Alpha = 12;end;
    if nargin<5,EnergyDistancePixelsize(1)=25; ... 
                EnergyDistancePixelsize(2)=.1; ...
                EnergyDistancePixelsize(3)= 0.36e-6;end;

% PARAMETERS.
folder      = '/home/moosmann/data/xenopus/';
filestring  = [folder 'xenopus000.edf'];
energy      = EnergyDistancePixelsize(1); % in keV
lambda      = EnergyConverter(energy); % in metres
distance    = EnergyDistancePixelsize(2); % in metres
pixelsize   = EnergyDistancePixelsize(3); % in metres
paddingProp = 2;
MethodOrPadValueProp = 'replicate';
paddingReco = 1;
PrintDomains= 0;
PadValue    = 0;
counts_min  = 0;
counts_max  = 0;
delta       = 3.88*10^-7;
k           = 2*pi/lambda;
BinaryFilter= 0.01;

% PROGRAM.
phase = k*delta*edfread(filestring)';
phase = phase(500+(1:512),1500+(1:512));
[dimx dimy] = size(phase);
% Blurring.
if Blurring(1)>0,
    hsizex = Blurring(1);
    hsizey = Blurring(2);
    sigma  = Blurring(3);
    phase  = imfilter(phase,fspecial('gaussian',[hsizex hsizey],sigma),MethodOrPadValueProp);
end;
% Propagation.
if 1,
int   = Propagation(phase,EnergyDistancePixelsize,paddingProp,MethodOrPadValueProp,PrintDomains);
phase = phase - mean(phase(:));
% Poisson noise.
if PhotonCountsForPoissonNoise>0,
    counts = PhotonCountsForPoissonNoise;
    if PhotonCountsForPoissonNoise*max(int(:))<2^16,
        int = imnoise(uint16(counts*int),'poisson');
        counts_min = min(int(:));
        counts_max = max(int(:));
        int = 1/counts*double(int);
    else
        fprintf(1,'Caution: Counts exceed 2^16\n');
    end;
end;
intmin = min(int(:));
intmax = max(int(:));
fprintf(1,['Energy: %2ukeV, Distance: %4.2fm, ' ...
           'Counts: %5u [%5u,%5u], Intensity: [%7.5f,%7.5f] (%.5f)\n'], ... 
        energy,distance,PhotonCountsForPoissonNoise, ... 
        counts_min,counts_max,intmin,intmax,intmax-intmin);
% Phase retrieval: CTF and nonlinear TIE.
if 1,
phiCTF = RecoCTF(int,Alpha,EnergyDistancePixelsize,BinaryFilter,paddingReco,PadValue);
phiBRO = Reco(int,Alpha,lambda,distance,pixelsize,1,0,0,1);
end;
end;

% Create proper namestring and print/save images: Exact phase. Intensity. Retrieved phase.
namestring  = sprintf('_E%02u_z%04.2f_c%05u',energy,distance,PhotonCountsForPoissonNoise);
namestring  = regexprep(namestring,'\.','p');
%mexVolWrite(sprintf('%sphase_int_FTint_CTF_BRO%s.vol',folder,namestring),cat(3,phase,int,intFT,phiCTF,phiBRO),'float32');


