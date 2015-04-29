% Script to test discrete rotation vs continous rotation. Be careful when
% Poisson noise is included. The blurred sinogram is acquired from
% averaging over disrete sinogram. Noise levels need to be adjusted
% appropriately before averaging such that noise levels of discrete and
% continous sinogram are comparable. Again, be careful not to confuse
% detector noise levels and sample noise levels for simulating a sinogram.
%
% Written by Julian Moosmann, last version: 2013-12-05, modified:
% 2014-11-12

%% Parameter
% Number of volume pixels
N = 1024;
% Number of projections
M = 2100;
% Oversampling factor
overSamp = 2;
NoiseLevel = 10000;
%
recoType = 'FBP_CUDA';'SART_CUDA';'SIRT_CUDA';
iterations = 100;
OutputPath = '/mnt/tomoraid-LSDF/users/moosmann/phd/figures/ContRotSim';
if NoiseLevel > 0 
    OutputPath = sprintf('%s/noiseLevel%05u',OutputPath,NoiseLevel);
end
%'~/Pictures/Rotation_flyScan_vs_stepwise';
CheckAndMakePath(OutputPath);
% ROI
x = N-256+1:N;
y = round(0.43*N)+(1:256);

%% Create Shepp-Logan phantom
fprintf('\n Create phantom (oversampled and downsampled)\n')
dataOver = phantom(N*overSamp*2);
dataOver = Binning(dataOver);
dataOver = dataOver/size(dataOver,1);

%% Create sinograms
% Oversample sino
fprintf('\n Create oversampled and downsampled sinogram and data\n')
[sinoOver, angles] = astra_make_sino(dataOver, N*overSamp, M);
sinoOver = sinoOver/overSamp;
sinoDown = zeros(M,N);
dataDown = zeros(N,N);
for nn = 1:overSamp
    sinoDown = sinoDown + sinoOver(:,nn:overSamp:end);
    dataDown = dataDown + dataOver(nn:overSamp:end,nn:overSamp:end);
end
sinoDown = sinoDown/overSamp;
dataDown = dataDown/overSamp;

% Non-oversampled data and sino
%data0 = phantom(N)/size(data0,1);
%sino0 = astra_make_sino(data0, N, M);
% domain(data0), domain(sino0),

% Print domains
fprintf('\n')
domain(dataOver),domain(dataDown),fprintf('\n')
domain(sinoOver),domain(sinoDown),fprintf('\n')

for Nred = 100;%[700 300 100]
    %% Create discrete and continous rotation sinogram
    projInc = round(M/Nred); % use odd number, the ceil(projInc)/2 is centered within the angular range to be blurred over
    projIndFirst = ceil(projInc/2);
    % discete sino
    sinod = sinoDown(projIndFirst+1:projInc:end,:);
    % continous blurred sino
    sinoc = zeros(size(sinod));
    for nn = 1:projInc
        sinoc = sinoc + sinoDown(nn:projInc:end,:);
    end
    sinoc = sinoc/projInc;
    anglesd = angles(projIndFirst:projInc:end);
    
    %% Noise
    if NoiseLevel > 0
        sinoc0 = sinoc;
        sinod0 = sinoc;
        addNoise = @(im) mean(im(:))/NoiseLevel*double(imnoise(uint16(NoiseLevel/mean(im(:))*(im)),'poisson'));
        sinoc = addNoise(sinoc0);
        sinod = addNoise(sinod0);
    end
        
    %% Reconstruction
    
    % discrete
    fprintf('\n Discrete reco: %s\n',recoType)
    recd = astra_make_reco(sinod,anglesd,recoType,iterations);
    
    % continous superpostions of single recos
    fprintf('\n Continous reco using superposition: %s\n',recoType)
    recStack = zeros(N,N,projInc);
    for nn = 1:projInc
        recStack(:,:,nn) = astra_make_reco(sinoc,angles(nn:projInc:end),recoType,iterations);
    end
    recd2 = recStack(:,:,projIndFirst);
    recc = sum(recStack,3)/projInc;
    
    fprintf('\n Continous reco from duplicated sino: %s\n',recoType)
    sinoc2 = zeros(size(sinoDown));
    % continous one step reco
    for nn = 1:projInc
        sinoc2(nn:projInc:end,:) = sinoc;
    end
    recc2 = astra_make_reco(sinoc2,angles,recoType,iterations);
    
    %% Save images
    SubOutputPath = sprintf('%s/projInc%02u_numProj%04uof%04u',OutputPath,projInc,size(sinod,1),M);
    CheckAndMakePath(SubOutputPath);
    hroi = @(im) im(x,y);
    % Reco, full
    h = @(im,imInd,string) WriteImage(sprintf('%s/%02u_%s_sizeFull_rangeFull_proj%04u_%s',SubOutputPath,imInd,recoType(1:4),size(sinod,1),string),im,'png');
    h(recd,  111, 'sinoDisc_recoSingle')
    h(recd2, 112, 'sinoCont_recoSingle')
    h(recc,  113, 'sinoCont_recoSingleSuperpos')
    %h(recc2, 114, 'sinoContDupl_recoSingle')
    % Reco, full, adjusted dyn range
    h = @(im,imInd,string) WriteImage(sprintf('%s/%02u_%s_sizeFull_range01_proj%04u_%s',SubOutputPath,imInd,recoType(1:4),size(sinod,1),string),SetDynRange(im,[0 1]),'png');
    h(recd,  121, 'sinoDisc_recoSingle')
    h(recd2, 122, 'sinoCont_recoSingle')
    h(recc,  123, 'sinoCont_recoSingleSuperpos')
    %h(recc2, 124, 'sinoContDupl_recoSingle')    
    % Reco, ROI
    h = @(im,imInd,string) WriteImage(sprintf('%s/%02u_%s_sizeROI_rangeFull_proj%04u_%s',SubOutputPath,imInd,recoType(1:4),size(sinod,1),string),hroi(im),'png');
    h(recd,  131, 'sinoDisc_recoSingle')
    h(recd2, 132, 'sinoCont_recoSingle')
    h(recc,  133, 'sinoCont_recoSingleSuperpos')
    %h(recc2, 134, 'sinoContDupl_recoSingle')
    % Difference Maps
    herr = @(im) abs(im - dataDown);    
    dynRange = [0 0.9*min( [ max2(herr(recd)) max2(herr(recd2)) max2(herr(recc)) max2(herr(recc2)) ] ) ];
    h = @(im,imInd,string) WriteImage(sprintf('%s/%02u_%s_err_sizeFull__rangAdj__proj%04u_%s',SubOutputPath,imInd,recoType(1:4),size(sinod,1),string),SetDynRange(herr(im),dynRange),'png');
    h(recd,  211, 'sinoDisc_recoSingle')
    h(recd2, 212, 'sinoCont_recoSingle')
    h(recc,  213, 'sinoCont_recoSingleSuperpos')
    %h(recc2, 214, 'sinoContDupl_recoSingle')
    % Difference Maps, ROI    
    herr = @(im) abs(hroi(im) - hroi(dataDown));
    dynRange = [0 0.9*min( [ max2(herr(recd)) max2(herr(recd2)) max2(herr(recc)) max2(herr(recc2)) ] ) ];
    h = @(im,imInd,string) WriteImage(sprintf('%s/%02u_%s_err_sizeROI__rangAdj__proj%04u_%s',SubOutputPath,imInd,recoType(1:4),size(sinod,1),string),SetDynRange(herr(im),dynRange),'png');
    h(recd,  221, 'sinoDisc_recoSingle')
    h(recd2, 222, 'sinoCont_recoSingle')
    h(recc,  223, 'sinoCont_recoSingleSuperpos')
    %h(recc2, 224, 'sinoContDupl_recoSingle')

end