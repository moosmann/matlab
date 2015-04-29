% Similiar script as 'ContRotRecoSim', but now looped over the projection increment to make a movie
%
% Script to test discrete rotation vs continous rotation. Be careful when
% Poisson noise is included. The blurred sinogram is acquired from
% averaging over disrete sinogram. Noise levels need to be adjusted
% appropriately before averaging such that noise levels of discrete and
% continous sinogram are comparable. Again, be careful not to confuse
% detector noise levels and sample noise levels for simulating a sinogram.
%
% Written by Julian Moosmann, last version: 2013-12-10

%% Parameter
% Number of volume pixels
N = 1024;
% Number of projections
M = 2100;
% Oversampling factor
overSamp = 2;
%
recoType = 'FBP_CUDA';'SART_CUDA';'SIRT_CUDA';
iterations = 100;
OutputPath = '~/Pictures/Rotation_flyScan_vs_stepwise';
CheckAndMakePath(OutputPath);
% ROI
x = N-256+1:N;
y = round(0.43*N)+(1:256);

%% Create Shepp-Logan phantom
fprintf('\n Create phantom (oversampled and downsampled)\n')
dataOver = Binning(phantom(N*overSamp*2));

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
data0 = phantom(N);
sino0 = astra_make_sino(data0, N, M);

% Loop
mmMax = 210;
recStack = zeros(N,N,mmMax);
for mm = 1:mmMax
    %% Create discrete and continous rotation sinogram
    projInc = mm; % use odd number, the ceil(projInc)/2 is centered within the angular range to be blurred over
    projIndFirst = ceil(projInc/2);
    % discete sino
    
    % continous blurred sino
    a = mm:mm:M;
    projEnd = a(end);        
    anglesd = angles(projIndFirst:projInc:projEnd);
    
    %% Reconstruction    
    % discrete
    
    recStack(:,:,mm) = astra_make_reco(sinoDown(projIndFirst:projInc:projEnd,:),anglesd(:),recoType,iterations);
    
end
outPath = '/ufs/moosmann/Pictures/Rotation_flyScan_vs_stepwise/';
WriteVol(recStack,[outPath 'RecoSequenceofIncreasingProjInc_discrete'])
WriteVol(recStack(x,y,:),[outPath 'RecoSequenceOfIncreasingProjInc_discrete'])