clear all

%% Parameter 
fprintf('\n\n--------------------------------------------------------------\n')
switch 'xeno4cell'
    case 'xeno4cell'
        intPath = '/mnt/tomoraid-LSDF/users/moosmann/CWI_DATA/ESRF_MI1079_ID19_July2011_inlineTomo/int/Xenopus_4cell_20keV/int_tif';
        intExp = imread(sprintf('%s/int_0001.tif',intPath));
        EnergyDistancePixelsize = [20 0.945 .75e-6];
    case 'xeno27'
        intPath = '/mnt/tomoraid-LSDF/users/moosmann/CWI_DATA/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/int/Xenopus_inVivo/Jul29_15-10_urea_stage27p0_30p0keV_0700mm_15ms_0500proj_scantime20s_deadtime8min/tomo01/proj_edf';
        intExp = pmedfread(sprintf('%s/int_0250.edf',intPath))';
        intExp = intExp(401:end-400,401:end-400);
        EnergyDistancePixelsize = [30 0.700 1.1e-6];
    case 'lena'
        pha0 = 0.1*normat(double(imread('~/lena.tif')));
        pha0 = -padarray(pha0,[256 256],0);
        pha0 = FilterBlur(pha0,[3 3], 1);
        EnergyDistancePixelsize = [30 0.500 1e-6];
        intExp = Propagation2(pha0,1*-10^-2.5*pha0,EnergyDistancePixelsize,'',1);        
end
outputPrecision = 'double';
Padding = '';
printInfo = 0;
bft = 0.1;
% Padding
%intExp = padarray(intExp,[size(intExp)/2],'symmetric','both');
imSize = size(intExp);
% Disc of interest
[x, y] = meshgrid(linspace(-1,1,imSize(2)),linspace(-1,1,imSize(1)));
doi = x.^2 + y.^2 < 0.85;
intExpMean = mean(intExp(:));

domain(intExp)

%% Phase retrieval
intfft2 = fft2(intExp-1);
phase = @(epsilon) ifft2( PhaseFilterDual('qp',size(intExp),EnergyDistancePixelsize,epsilon,bft,outputPrecision).* intfft2);
%leps=1:0.2:10;for nn = numel(leps):-1:1,pha(:,:,nn) = phase(10^-leps(nn));end

%% Simulation of intensity
roi = [284  1437  105  158];
% Attenuation
BctfMean = -mean(intExp(:)-1)/2;
fprintf(' Dual CTF: mean(B): %g, exp(-2*mean(B): %g\n',BctfMean,exp(-2*BctfMean))
BMeanExp = -log(intExpMean)/2;
fprintf(' mean(intExp): %g, -ln(mean(intExp))/2: %g\n',intExpMean,BMeanExp)
M = 10;
% Duality
N = 10;
lepsMin = 0.0;
lepsMax = 5;
logepsilon = -(lepsMin + (lepsMax-lepsMin).*(0:N)/N);
epsilon = 10.^logepsilon;
% Loop
fprintf('\n')

for nn = numel(epsilon):-1:1
    % Phase
    pha = phase(epsilon(nn));
    % set
    pha = pha -max(pha(:));
    phaStack(:,:,nn) = pha;
    % Simulate intensity
    %int = Propagation2(pha,BMean(mm)-epsilon(nn)*pha,EnergyDistancePixelsize,Padding,printInfo);
    int = Propagation2(pha,-epsilon(nn)*pha,EnergyDistancePixelsize,Padding,printInfo);
    intStack(:,:,nn) = int;
    %intStack(:,:,nn,mm) = int;
    % Measures
    intMean(nn) = mean( int(:) );    
    betaVar = 0.1;
    beta = intMean(nn)/intExpMean*( 1 + betaVar*(-M:2:M)/M );
    for mm = numel(beta):-1:1
        %intDiffNorm(nn,mm) = norm( int - intExp*exp(-2*BMean(mm)) ,1);
        epsilonGrid(nn,mm) = epsilon(nn);
        alpha(nn,mm) = log(beta(mm));
        intDiffNorm(nn,mm) = norm( (int - intExp*beta(mm)), 1);
        %h = @(x) fft2((x));
        %intDiffNorm(nn,mm) = norm( log(1+abs( h(int) - h(intExp*beta(mm)) )), 1);
        %intDiffROINorm(nn,mm) = norm( int(roi(1)+(1:roi(3)), roi(2)+(1:roi(4)), n) ,1);
        % Print info
        %domain(int,1,sprintf('nn:%2u eps:%10g',nn,epsilon(nn)));
        intMin = min(int(:));
        intMax = max(int(:));
        fprintf(' nn:%3u, mm:%3u, beta:%4.02g, eps:%6.3f, -log(eps):%5.3g  MEAN: Sim:%5.3g, Exp:%5.3g  MAX-MIN: int: %5.3g, eps*pha:%5.3g  intDiffNorm:%5.1f\n', ...
            nn,mm,beta(mm),epsilon(nn),-logepsilon(nn),intMean(nn),intExpMean*beta(mm),intMax-intMin,epsilon(nn)*(max(pha(:))-min(pha(:))),intDiffNorm(nn,mm))
    end
end
fprintf('\n')
%% Print and plot information
domain(intExp)
[~,intDiffNormMinPos]=min(intDiffNorm(:));
[nn, mm] = ind2sub(size(intDiffNorm),intDiffNormMinPos);
fprintf(' Minimum: nn: %u, mm: %u, beta: %g, -log10(epsilon): %g, epsilon: %g\n',nn,mm,beta(mm),-logepsilon(nn),epsilon(nn))
outputPath = '~/matlabWS/';
%save(sprintf('%s/DetOfEpsilon_%s.mat',outputPath,datestr(now,'yyyy-mm-dd_HH-MM-SS')));
ishow(intDiffNorm)
figure('Name','Contour over index: L1-Norm of difference of simulated and experimental intensity'),
contourf(flipud(intDiffNorm))
%figure('Name','Contour over epsilon & alpha: L1-Norm of difference of simulated and experimental intensity')
%contourf(flipud(epsilonGrid),flipud(alpha),flipud(intDiffNorm))
fprintf('--------------------------------------------------------------\n')
