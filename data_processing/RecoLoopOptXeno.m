function RecoLoopOptXeno(alpha,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)
% Script for phase retrieval of Xenopus frog embryo at a 4 cell stage. Data
% set was taken at ESRF in Juli 2011. Script loops over all folders in the
% given parent folder and starting with Xenopus.

%% Default parameters.
if nargin < 1
    alpha        = [2.5 2.5]; %Regularization parameter
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
    BinaryFilter = 0.30;%if >0: evaluate projected CTF (PCTF), where value is the threshold for the projection filter.
end
if nargin < 6
    DataSetsToProcess = [29 27];
end
if nargin < 7
    %[Energy Distance Pixelsize] in metre.
%     offset = 0.023;
%     EnergyDistancePixelsize{1} = [22 offset+0.100 0.75e-6];
%     EnergyDistancePixelsize{2} = [22 offset+0.100 0.75e-6];
%     EnergyDistancePixelsize{3} = [22 offset+0.400 0.75e-6];
%     EnergyDistancePixelsize{4} = [22 offset+0.400 0.75e-6];
end
if nargin < 8
    ParentPath = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo/';
end
if nargin < 9
    DataFolderNamePrefix = 'opt_xeno_';
end
% 0000000001111111111222222222233333333333444444
% 1234567890123456789012345678901234567890123456
% opt_xeno_stage11_tomo_22keV_10x_125mm_0p20sec_
DataFolderStruct = dir([ParentPath 'data/' DataFolderNamePrefix '*']);
DistanceOffset = 40; %mm
for nn = numel(DataFolderStruct):-1:1
    Energy    =  str2double(DataFolderStruct(nn).name(23:24));
    Distance  = (DistanceOffset + str2double(DataFolderStruct(nn).name(33:35)))/1000;
    Pixelsize = 15e-6/str2double(DataFolderStruct(nn).name(29:30));
    EnergyDistancePixelsize{nn} = [Energy Distance Pixelsize];
    %fprintf('%2u. %-50s: Energy = %2u Distance = %5.3f Pixelsize = %7.2g\n',nn,DataFolderStruct(nn).name,EnergyDistancePixelsize{nn})
end

%% Call reconstruction loop.
RecoLoop(alpha,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilter,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)