function [out phase int] = supernova()

%% Parameters
OutputDir = '/home/moosmann/supernova/';
EnergyDistancePixelsize = [30 0.62 1e-6];%[20 1 1e-6];
prefac(EnergyDistancePixelsize(1),EnergyDistancePixelsize(2),EnergyDistancePixelsize(3));
blurring=[8 8 2];
filestring='/home/moosmann/lena.tif';
PhaseMin = 0.001;
PhaseMax = 10+PhaseMin;
NumStep  = 50;
PhaseRange = PhaseMin + (PhaseMax-PhaseMin)*(0:NumStep-1)/(NumStep-1);
PhaseRange = PhaseRange(1)+normat(PhaseRange.^3)*(PhaseRange(end)-PhaseRange(1));
 PhaseRange = [0.01 1 3 6];
% NumStep = numel(PhaseRange);
fprintf('Phase range: %s\n',mat2str(PhaseRange))
NoiseLevel = 10000;
%% Create phase map
% phase = zeros(1024);
% phase(257:768,257:768) = normat(double(imread(filestring)));
 phase = 0.0 + 1*normat(double(imread(filestring)));
 [dimx dimy] = size(phase);
phase = padarray(phase,[dimx/2 dimy/2],'symmetric','both');
% Gaussian blur the image to smooth the edges
if blurring(1) > 0
    hsizex = blurring(1);
    hsizey = blurring(2);
    sigma  = blurring(3);
    phase  = imfilter(phase,fspecial('gaussian',[hsizex hsizey],sigma),'symmetric');
end;
phase = normat(phase);
% Add noise
%% Propagate the phase map
for nn = numel(PhaseRange):-1:1
    int(:,:,nn) = Propagation(PhaseRange(nn)*phase,EnergyDistancePixelsize,1,'symmetric',0);
    %int = 1/NoiseLevel*double(imnoise(uint16(NoiseLevel*int),'poisson'));
    %int = int - mean(int(:));
    intft = (fftshift(fft2(int(:,:,1))));
    intfta = abs(intft);
    fprintf('phase max: %7.3g. ',PhaseRange(nn)),domain(intfta);
    out(:,:,nn) = log(0.0001+intfta);
    %domain(out(:,:,nn))   
end
% Normalize stack
outMin = min(out(:));
outMax = max(out(:));
out = (out-outMin)./(outMax-outMin);
%% Save images
OutputDirTif = [OutputDir 'tif/'];
OutputDirEdf = [OutputDir 'edf/'];
if ~exist(OutputDirTif,'dir')
    mkdir(OutputDirTif)
end
% if ~exist(OutputDirEdf,'dir')
%     mkdir(OutputDirEdf)
% end
% for nn = numel(PhaseRange):-1:1
%     imwrite(out(:,:,nn),sprintf('%slogabsftint_%02u.tif',OutputDirTif,nn))
%     edfwrite(sprintf('%slogabsftint_%02u.edf',OutputDirEdf,nn),out(:,:,nn),'float32');
% end