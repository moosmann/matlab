function RecoLoopXenopus4cell(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)
% Script for phase retrieval of Xenopus frog embryo at a 4 cell stage. Data
% set was taken at ESRF in Juli 2011. For details see function 'RecoLoop'.

%% Default parameters.
if nargin < 1
    %Regularization parameter
    alphaCTF_alphaTIE = 2.5;
end
if nargin < 2
    %Evaluate leading-order (lo) transport-of-intensity (TIE) expression TIElo
    evalTIElo    = 0;
end
if nargin < 3
    %Evaluate perturbatively next-to-leading-order (nlo) TIE expression: TIEpnlo
    evalTIEpnlo  = 0;
end
if nargin < 4
    %Evaluate contrast-transfer-function (CTF) expression.
    evalCTF      = 0;
end
if nargin < 5
    %if >0: evaluate projected CTF (PCTF), where value is the threshold for the projection filter.
    BinaryFilter = 7.1;
end
if nargin < 6
    DataSetsToProcess = 1;
end
if nargin < 7
    %[Energy Distance Pixelsize] in metre.
    EnergyDistancePixelsize{1} = [20 0.945 0.745e-6];
    EnergyDistancePixelsize{2} = [40 0.945 0.745e-6];
end
if nargin < 8
    %ParentPath = '/mnt/tomoraid3/user/moosmann/Xenopus_4cell/';
    ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell';
end
if nargin < 9
    DataFolderNamePrefix = 'Xenopus';
end

%% Call reconstruction loop.
RecoLoop(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)
