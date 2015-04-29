function RecoLoopZebrafish(alpha,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)
% Script for phase retrieval of Xenopus frog embryo at a 4 cell stage. Data set was taken at
% ESRF in Juli 2011. Script loops over all folders in the given parent
% folder and starting with Xenopus.

%% Default parameters.
if nargin < 1
    alpha        = 2.5; %Regularization parameter
end
if nargin < 2
    evalTIElo    = 1;%Evaluate leading-order (lo) transport-of-intensity (TIE) expression TIElo
end
if nargin < 3
    evalTIEpnlo  = 1;%Evaluate perturbatively next-to-leading-order (nlo) TIE expression: TIEpnlo
end
if nargin < 4
    evalCTF      = 0;%Evaluate contrast-transfer-function (CTF) expression.
end
if nargin < 5
    BinaryFilter = 0.1;%if >0: evaluate projected CTF (PCTF), where value is the threshold for the projection filter.
end
if nargin < 6
    DataSetsToProcess = 0;
end
if nargin < 7
    %[Energy Distance Pixelsize] in metre.
    EnergyDistancePixelsize = [20 0.945 0.745e-6];
end
if nargin < 8
    %ParentPath = '/mnt/tomoraid3/user/moosmann/Zebrafish_MI1079-ESRF-ID19-July2011_inlineTomo/';
    ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/ZebraFish';
end
if nargin < 9
    DataFolderNamePrefix = 'ZebraFish';
end

%% Call reconstruction loop.
RecoLoopGPU(alpha,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)