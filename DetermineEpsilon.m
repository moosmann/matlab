clear all

%% Parameter 
fprintf('\n\n--------------------------------------------------------------\n')
switch 'xeno4cell'
    case 'xeno4cell'
        intPath = '/mnt/tomoraid-LSDF/users/moosmann/CWI_DATA/ESRF_MI1079_ID19_July2011_inlineTomo/int/Xenopus_4cell_20keV/int_tif';
        intExp = imread(sprintf('%s/int_0001.tif',intPath));
        EnergyDistancePixelsize = [20 0.945 .75e-6];
    case 'xeno27'
        intPath = '/mnt/tomoraid-LSDF/users/moosmann/CWI_DATA/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/int/Xenopus_inVivo/Jul29_15-10_urea_stage27p0_30p0keV_0700mm_15ms_0500proj_scantime20s_deadtime8min/tomo01';
        intExp = pmedfread(sprintf('%s/int_0250.edf',intPath))';
        intExp = intExp(401:end-400,401:end-400);
        EnergyDistancePixelsize = [30 0.700 1.1e-6];
    case 'lena'
        pha0 = 1*normat(double(imread('~/lena.tif')));
        pha0 = -padarray(pha0,[256 256],0);
        pha0 = FilterBlur(pha0,[3 3], 1);
        EnergyDistancePixelsize = [30 8*0.945 1.1e-6];
        intExp = Propagation2(pha0,-10^-3*pha0,EnergyDistancePixelsize,'',1);
end
%% Parameter
outputPrecision = 'double';
Padding = '';
printInfo = 0;
bft = 0.01;
imSize = size(intExp);
% Disc of interest
[x, y] = meshgrid(linspace(-1,1,imSize(2)),linspace(-1,1,imSize(1)));
doi = x.^2 + y.^2 < 0.85;
intExpMean = mean(intExp(:));
domain(intExp)
%intExp = intExp/intExpMean;
%intExpMean = mean(intExp);
%domain(intExp)

%% Phase retrieval
gfft2 = fft2(intExp-1);
phase = @(epsilon) ifft2( PhaseFilterDual('qp',size(intExp),EnergyDistancePixelsize,epsilon,bft,outputPrecision).* gfft2);

%% Simulation of intensity
% Attenuation
BMeanCTF = -mean(intExp(:)-1)/2;
BMeanExp = -log(intExpMean)/2;
fprintf(' Dual CTF: mean(B)=-mean(intExp-1)/2: %g, exp(-2*mean(B)): %g\n',BMeanCTF,exp(-2*BMeanCTF))
fprintf(' mean(intExp): %g, -ln(mean(intExp))/2: %g\n',intExpMean,BMeanExp)
% Duality
N = 20;
lepsMin = 1;
lepsMax = 4;
logepsilon = -(lepsMin + (lepsMax-lepsMin).*(0:N)/N);
epsilon = 10.^logepsilon;
% beta
M = 10;
beta = 1 + 0.1*(-M:2:M)/M; 
% Loop
fprintf('\n')
for nn = numel(epsilon):-1:1
    %% Compute phase
    pha = phase(epsilon(nn));
    % set
    pha = pha -max(pha(:));    
    phaStack(:,:,nn) = pha;
    %% Simulate intensity
    int = Propagation2(pha,-epsilon(nn)*pha,EnergyDistancePixelsize,Padding,printInfo);    
    intStack(:,:,nn) = int;
    intMean(nn) = mean(int(:));    
    for mm = numel(beta):-1:1               
        %% Difference Norm
        p = 1;
        h = @(x) sum( abs( x(:) ).^p )^1/p;        
        intDiffNorm(nn,mm) = h( beta(mm)*int/intMean(nn) - intExp/intExpMean );
        
        % Print info
        intMin = min(int(:));
        intMax = max(int(:));
        fprintf(' nn:%3u, mm:%3u, beta:%4.02g, eps:%6.3f, -log(eps):%5.3g. INT: Mean:%4.2g, Max-Min:%4.2g, Min:%4.2g, Max:%4.2g, intDiffNorm:%5.1f\n', ...
            nn,mm,beta(mm),epsilon(nn),-logepsilon(nn),intMean(nn),intMax-intMin,intMin,intMax,intDiffNorm(nn,mm))        
    end
end
fprintf('\n')

%% Find minimum of Difference Norm
domain(intExp)
[~, intDiffNormMinPos] = min(intDiffNorm(:));
[nn, mm] = ind2sub(size(intDiffNorm),intDiffNormMinPos);
fprintf(' Minimum: nn: %u, mm: %u, beta: %g, -log10(epsilon): %g, epsilon: %g, intDiffNorm:%5.1f\n',nn,mm,beta(mm),-logepsilon(nn),epsilon(nn),intDiffNorm(nn,mm))

%% Print and plot information
outputPath = '~/matlabWS/';
%save(sprintf('%s/DetOfEpsilon_%s.mat',outputPath,datestr(now,'yyyy-mm-dd_HH-MM-SS')));
% PLOTS
%ishow(intDiffNorm)
figure('Name','Contour over index: L1-Norm of difference of simulated and experimental intensity'),
contourf(flipud(intDiffNorm))


%figure('Name','epsilon  VS  intDiffMean'),plot(epsilon(:),intDiffNorm(:),'.')
%figure('Name','log10(epsilon)  VS  intDiffMean'),plot(-logepsilon(:),intDiffNorm(:),'.')

fprintf('--------------------------------------------------------------\n')
