% Load Lena test pattern, pad it, and optionally blur it with a
%  Gaussian. Compute intensity and optionally noise it.
clear all
%% Parameters
blurring = 0;%[12 12 1];%[16 16 16];
Energy = 20;%keV
Distance = 0.5;%m
Pixelsize = 1e-6;%m
NoiseLevel = 000;%20000;
filestring='/home/moosmann/lena.tif';
OutputPath = '/home/moosmann/simulation/';
maxPhaShifts = [0.01 0.1 1 6];%[0.01 0.1 1 3 6];
Padding = 1;
PadMethodOrValue = 'symmetric';
alpha = 10;
PaddingFacMeth = {2 'symmetric'};
%% Read lena test pattern to create a phase map.
% phase = zeros(1024);
% phase(257:768,257:768) = normat(double(imread(filestring)));
phase = normat(double(imread(filestring)));
phase = padarray(phase,size(phase)/2,'symmetric','both');
%% Propagate the image.
for ii = numel(maxPhaShifts):-1:1
    int(:,:,ii) = Propagation(maxPhaShifts(ii)*phase,[Energy Distance Pixelsize],Padding,PadMethodOrValue,0);
    domain(int(:,:,ii),1,sprintf('intensity for phi_max = %g',maxPhaShifts(ii)))
end
%% Gaussian blur intensities
if blurring(1) > 0
    for ii = size(int,3):-1:1
        hsizex = blurring(1);hsizey = blurring(2);sigma  = blurring(3);
        int(:,:,ii)  = imfilter(int(:,:,ii),fspecial('gaussian',[hsizex hsizey],sigma),'symmetric','same','conv');
    end
end;
NumIm = size(int,3);
%% FT of intensity
for ii = NumIm:-1:1
    m                = fft2(int(:,:,ii));
    intFT(:,:,ii)    = m;
    m(1,1) = 0;
    intFTabsLOG(:,:,ii) = log(1+abs(fftshift(m)));
end
%% Crop images
phase = phase(257:768,257:768);
int = int(257:768,257:768,:);
%% phase retrieval
for ii = NumIm:-1:1
    out(ii) = Reco(int(:,:,ii),alpha,[Energy Distance Pixelsize],[1 0.1 1 0],0,PaddingFacMeth);
end
%% Show movies
% NormGlobally = 0;
% nimplay(int,NormGlobally)
% nimplay(intFTabsLOG,NormGlobally)
% nimplaystruct(out,'tieLO',NormGlobally)
% nimplaystruct(out,'ctf',NormGlobally)
% nimplaystruct(out,'ctfProjected',NormGlobally)
%% Save images
for ii = NumIm:-1:1
    edfwrite(sprintf('%sint_%04u_phimax%04g.edf',OutputPath,ii,maxPhaShifts(ii)),int(:,:,ii)','float32');
    edfwrite(sprintf('%sintFTabsLOG_%04u_phimax%04g.edf',OutputPath,ii,maxPhaShifts(ii)),intFTabsLOG(:,:,ii)','float32');
    edfwrite(sprintf('%sphaseTIE_%04u_phimax%04g.edf',OutputPath,ii,maxPhaShifts(ii)),out(ii).tieLO','float32');
    edfwrite(sprintf('%sphaseCTF_%04u_phimax%04g.edf',OutputPath,ii,maxPhaShifts(ii)),out(ii).ctf','float32');
    edfwrite(sprintf('%sphasePCTF_%04u_phimax%04g.edf',OutputPath,ii,maxPhaShifts(ii)),out(ii).ctfProjected','float32');
end
edfwrite(sprintf('%sphase.edf',OutputPath),phase','float32');
%% Noise
if NoiseLevel > 0
    %% Add noise to the intensity.
    for ii = NumIm:-1:1
        intPN(:,:,ii) = 1/NoiseLevel*double(imnoise(uint16(NoiseLevel*int(:,:,ii)),'poisson'));
    end
    intPN = intPN(257:768,257:768,:);
    %% FT of intensity
    for ii = NumIm:-1:1
        m                  = fft2(intPN(:,:,ii));
        intPNFT(:,:,ii)    = m;
        m(1,1) = 0;
        intPNFTabs(:,:,ii) = abs(fftshift(m));
        outPN(ii) = Reco(intPN(:,:,ii),alpha,[Energy Distance Pixelsize],[1 0.1 1 0],0,PaddingFacMeth);
    end
    %% Show movies
    nimplay(intPN,NormGlobally)
    nimplay(log(1+intPNFTabs),NormGlobally)
    nimplaystruct(out,'tieLO',NormGlobally)
    nimplaystruct(out,'ctf',NormGlobally)
    nimplaystruct(out,'ctfProjected',NormGlobally)
end