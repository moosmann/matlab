function [phase,int,phiCTF,phiBRO,EnergyDistancePixelsize,intFT]=Container(PhotonCountsForPoissonNoise,compute_correction,alpha,EnergyDistancePixelsize)

% function [phase,int,phiCTF,phiBRO,EnergyDistancePixelsize,intFT]=Container(PhotonCountsForPoissonNoise,compute_correction,alpha,EnergyDistancePixelsize)

%% DEFAULTS.
if nargin<1,PhotonCountsForPoissonNoise=0;end;
if nargin<2,compute_correction = 0;end;
if nargin<3,alpha = 12;end;
if nargin<4,EnergyDistancePixelsize(1)=25; ...
        EnergyDistancePixelsize(2)=1; ...
        EnergyDistancePixelsize(3)= 0.36e-6;
end;
%% PARAMETERS.
folder      = '/home/moosmann/data/xenopus/';
filestring  = [folder 'container000.edf'];
energy      = EnergyDistancePixelsize(1); % in keV
lambda      = EnergyConverter(energy); % in metres
distance    = EnergyDistancePixelsize(2); % in metres
pixelsize   = EnergyDistancePixelsize(3); % in metres
paddingProp = 2;
MethodOrPadValueProp = 'replicate';
paddingReco = 1;
PrintDomains= 1;
PadValue    = 0;
counts_min  = 0;
counts_max  = 0;
delta       = 3.88*10^-7;
k           = 2*pi/lambda;
%% PROGRAM.
phase = k*delta*edfread(filestring)';
domain(phase,'imported image');
%[dimx dimy] = size(phase);
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
            fprintf(1,'Warning: Counts exceed 2^16\n');
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
        phiCTF = RecoCTF(int,alpha,lambda,distance,pixelsize,paddingReco,PadValue,compute_correction);
        phiBRO = Reco(int,alpha,lambda,distance,pixelsize,1,0,0,1);
        intFT  = laf(int);
    end;
end;

