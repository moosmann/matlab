function RecoLoopXenopusESRF(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)
% Script for phase retrieval. Script loops over all folders under the parent
% folder 'ParentPath' and starting with 'Xenopus'

%% Default parameters.
if nargin < 1
    alphaCTF_alphaTIE = 2.5;%[2.5 2.5]; %Regularization parameter
end
if nargin < 2
    evalTIElo = 0;%Evaluate leading-order (lo) transport-of-intensity (TIE) expression TIElo
end
if nargin < 3
    evalTIEpnlo = 0;%Evaluate perturbatively next-to-leading-order (nlo) TIE expression: TIEpnlo
end
if nargin < 4
    evalCTF = 0;%Evaluate contrast-transfer-function (CTF) expression.
end
if nargin < 5
    BinaryFilter = 0.1;%if >0: evaluate projected CTF (PCTF), where value is the threshold for the projection filter.
end
if nargin < 6
    DataSetsToProcess = 1;%[4 3 2 5 1];
end
if nargin < 7
    %[Energy Distance Pixelsize] in metre.
    EnergyDistancePixelsize = [20 0.94865 1.4e-6];
end
if nargin < 8
    ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_May2011_Xenopus/';
end
if nargin < 9
    DataFolderNamePrefix = 'Xenopus';
end

RecoLoop(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)
