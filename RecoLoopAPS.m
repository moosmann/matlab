function [sinotie sinopctf] = RecoLoopAPS(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,EnergyPixelsizeDistance,InputPath,FileNamePrefix,OutputPath,doSino,IndexArray,DataSetName)
% Script for phase retrieval of Xenopus frog embryos. Data set was taken at
% ESRF in May 2011. Script loops over all folders contained in the CURRENT 'parent
% folder' and starting with Xenopus. Runs on laptotp. Makes sinogram of
% central slice. hdf reading and writing not needed. In general consistency
% has to be checked. Called function RecoLoopThis will be changed.

%% Default parameters.
if nargin < 1
    alphaCTF_alphaTIE  = 2.5;
end
if nargin < 2
    evalTIElo    = 1;
end
if nargin < 3
    evalTIEpnlo  = 0;
end
if nargin < 4
    evalCTF      = 0;
end
if nargin < 5
    BinaryFilterThreshold = 0.1;
end
if nargin < 6
    EnergyPixelsizeDistance = [20 1.4e-6 0.94865];
end
if nargin < 7
    InputPath = pwd; % String: trailing seperator ('/') not needed. E.g.
end
if nargin < 8
    FileNamePrefix = 'int';
end
if nargin < 9
    OutputPath = [pwd '/phase/'];
end
if nargin < 10
    doSino = 1;
end
if nargin < 11
     IndexArray  = {[1  1],[1  1],[1024  2048]};
end
if nargin < 12
    DataSetName = '/entry1/data/data';
end
%% Call RecoLoopThis
[sinotie sinopctf] = RecoLoopThis(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,EnergyPixelsizeDistance,InputPath,FileNamePrefix,OutputPath,doSino,IndexArray,DataSetName);