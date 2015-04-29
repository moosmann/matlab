%Script for quick phase retrieval

%% Parameters for reading the data.
DataPath = '/data/visitor/mi1057/bm05';
DataType = 'edf';
FirstFileToRead = 1;
StepSize = 1;
NumSteps = 1;
%% Parameters for phase retrieval.
alpha = [2.5 5];
EnergyDistancePixelsize = [25 1.5e-6];
evalTIElo = 1;
evalTIEpnlo = 0;
evalCTF = 0;
BinaryFilterThreshold = 0.01;
%% Data preprocessing.
% Read data.
dat = readstack(DataPath,'Xenopus*0000',DataType,FirstFileToRead,StepSize,NumSteps);
NumOfImages = size(dat,3);
% Read refs.
ref = readstack(DataPath,'ref*_0000',DataType,FirstFileToRead,StepSize,NumSteps);
if size(ref,3) > 1
    ref = median(ref,3);
end
% Optionally read dark fields.
DarkAr = dir([DataPath '/*dark*.edf']);
% Flat-(and dark-)field correction.
if numel(DarkAr) == 1
    dark = readstack(DataPath,DarkAr(1).name,DataType,FirstFileToRead,StepSize,NumSteps);
    im = (im - repmat(dark,[1 1 NumOfImages]))./(repmat(ref,[1 1 NumOfImages]) - repmat(dark,[1 1 NumOfImages]));
else
    im = (im)./(ref);
end
%% Phase retrieval.
for ii = size(im,3):-1:1
    out(ii) = Reco(im(ii),alpha,EnergyDistancePixelsize,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,Padding_FactorAndValue);
end


