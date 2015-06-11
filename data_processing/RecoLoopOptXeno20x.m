function RecoLoopOptXeno20x(alpha,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)
% Script for phase retrieval of Xenopus frog embryo at a 4 cell stage. Data
% set was taken at ESRF in Juli 2011. Script loops over all folders in the
% given parent folder and starting with Xenopus.

%% Default parameters.
if nargin < 1
    alpha        = [5 2.5]; %Regularization parameter
end
if nargin < 2
    evalTIElo    = 0;%Evaluate leading-order (lo) transport-of-intensity (TIE) expression TIElo
end
if nargin < 3
    evalTIEpnlo  = 0;%Evaluate perturbatively next-to-leading-order (nlo) TIE expression: TIEpnlo
end
if nargin < 4
    evalCTF      = 0;%Evaluate contrast-transfer-function (CTF) expression.
end
if nargin < 5
    BinaryFilter = 0.01;%if >0: evaluate projected CTF (PCTF), where value is the threshold for the projection filter.
end
if nargin < 6
    DataSetsToProcess = [1 3 2 4];
end
if nargin < 7
    %[Energy Distance Pixelsize] in metre.
    offset = 0.023;
    EnergyDistancePixelsize{1} = [22 offset+0.100 0.75e-6];
    EnergyDistancePixelsize{2} = [22 offset+0.100 0.75e-6];
    EnergyDistancePixelsize{3} = [22 offset+0.400 0.75e-6];
    EnergyDistancePixelsize{4} = [22 offset+0.400 0.75e-6];
end
if nargin < 8
    ParentPath = '/mntdirect/_data_visitor/mi1057/bm05/pc/opt_xeno_stage12_tomo';
end
if nargin < 9
    DataFolderNamePrefix = 'xeno';
end

%% Call reconstruction loop.
RecoLoop(alpha,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)