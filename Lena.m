function [phase, int, int0] = Lena(MaxPhaseShift,blurring,EnergyDistancePixelsize,AverageCounts,filestring)

% Load Lena test pattern, pad it, and optionally blur it with a
%  Gaussian. Compute intensity and optionally noise it.
  
if nargin < 1
    MaxPhaseShift = 1;
end
if nargin < 2
    blurring=[16 16 16];
end;
if nargin < 3
    EnergyDistancePixelsize = [20 .7 1.1e-6];
end
if nargin < 4
    AverageCounts = 20000;
end
if nargin < 5
    filestring='~/lena.tif';
end;

%% Read lena test pattern to create a phase map.
phase = zeros(1024);
phase(257:768,257:768) = MaxPhaseShift*normat(double(imread(filestring)));
%% Gaussian blur the image to smooth the edges.
if blurring(1) > 0
    hsizex = blurring(1);
    hsizey = blurring(2);
    sigma  = blurring(3);
    phase  = imfilter(phase,fspecial('gaussian',[hsizex hsizey],sigma));
end;
phase = MaxPhaseShift*normat(phase);
%% Propagate the image.
int = Propagation(phase,EnergyDistancePixelsize,2,'symmetric',0);
if nargout == 3
    int0 = int;
end
%phase = phase - mean(phase(:));
%% Add noise to the intensity.
if AverageCounts > 0
    int = 1/AverageCounts*double(imnoise(uint16(AverageCounts*int),'poisson'));
end;
