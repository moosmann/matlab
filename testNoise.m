
%% Read lena test pattern to create a phase map.
filestring='/home/moosmann/lena.tif';
phase = zeros(1024);
phase(257:768,257:768) = normat(double(imread(filestring)));
NoiseLevel = [ 1 10 100 200 500 1000 2000 5000 10000 20000 50000];
%% Add noise
for ii = numel(NoiseLevel):-1:1
        phasePN(:,:,ii) = 1/NoiseLevel(ii)*double(imnoise(uint16(NoiseLevel(ii)*phase),'poisson'));
        domain(phasePN(:,:,ii))
        %intGN = 1/NoiseLevel(ii)*double(imnoise(uint16(NoiseLevel(ii)*phase),'poisson'));
end
%% 
phaseDiff = phasePN - repmat(phase,[1 1 numel(NoiseLevel)]);
%% Plots
nimplay(phasePN)
nimplay(phaseDiff)

