function RecoLoop2(alpha,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,DataSetsToProcess,edp,DataPath,SubDataPathPrefix)
% Script for phase retrieval of Xenopus frog embryos. Data set was taken at
% ESRF in May 2011. Script loops over all folders in the given parent
% folder and starting with Xenopus.

%% Default parameters.
if nargin < 1
    alpha        = 2.5; %Regularization parameter
end
if nargin < 2
    evalTIElo    = 1;%Evaluate leading-order (lo) transport-of-intensity (TIE) expression TIElo
end
if nargin < 3
    evalTIEpnlo  = 0;%Evaluate perturbatively next-to-leading-order (nlo) TIE expression: TIEpnlo
end
if nargin < 4
    evalCTF      = 0;%Evaluate contrast-transfer-function (CTF) expression.
end
if nargin < 5
    BinaryFilterThreshold = 0.01;%if >0: evaluate projected CTF (PCTF), where value is the threshold for the projection filter.
end
if nargin < 6
    DataSetsToProcess = 0;
end
if nargin < 7
    edp          = [20 0.94865 1.4e-6];%[Energy Distance Pixelsize] in metre.1
end
if nargin < 8
    DataPath = '/mnt/tomoraid3/user/moosmann/Xenopus_ESRF_May2011/';
end
if nargin < 9
    SubDataPathPrefix = 0;
end

%% Read folder names of all data sets.
if SubDataPathPrefix == 0
    SubDataPathPrefix     = DataPath(1:end-1);
    [~,SubDataPathPrefix] = fileparts(SubDataPathPrefix);
    SubDataPathPrefix     = SubDataPathPrefix(1:4);
end
DataFolderNames = dir([DataPath 'data/' SubDataPathPrefix '*']);
inputParentPath  = [DataPath 'int/'];
outputPath = [DataPath 'phase/'];
%% Loop over different data sets (frog embryo stages).
tic
totalcounter = 0;
if DataSetsToProcess == 0
    DataSetsToProcess = 1:numel(DataFolderNames);
end
for kk = DataSetsToProcess
    % Create folders for phase images.
    if evalTIElo
        outputFolderTIElo = sprintf('%s/TIElo_alpha%3.2f',DataFolderNames(kk).name,alpha);
        outputFolderTIElo = regexprep(outputFolderTIElo,'\.','p');
        mkdir(outputPath,outputFolderTIElo);
        outputPathTIElo   = [outputPath outputFolderTIElo '/'];
    end
    if evalTIEpnlo
        outputFolderTIEpnlo = sprintf('%s/TIEpnlo_alpha%3.2f',DataFolderNames(kk).name,alpha);
        outputFolderTIEpnlo = regexprep(outputFolderTIEpnlo,'\.','p');
        mkdir(outputPath,outputFolderTIEpnlo);
        outputPathTIEpnlo   = [outputPath outputFolderTIEpnlo '/'];
    end
    if evalCTF
        outputFolderCTF = sprintf('%s/CTF_alpha%3.2f',DataFolderNames(kk).name,alpha);
        outputFolderCTF = regexprep(outputFolderCTF,'\.','p');
        mkdir(outputPath,outputFolderCTF);
        outputPathCTF   = [outputPath outputFolderCTF '/'];
    end
    if BinaryFilterThreshold(1) > 0
        if size(BinaryFilterThreshold,2) == 1
            outputFolderCTFproj = sprintf('%s/CTFproj_alpha%3.2f_binFilt%3.4f',DataFolderNames(kk).name,alpha,BinaryFilterThreshold(1));
        else
            outputFolderCTFproj = sprintf('%s/CTFproj_alpha%3.2f_binFilt%3.4fblurW%02uS%01u',DataFolderNames(kk).name,alpha,BinaryFilterThreshold);
        end
        outputFolderCTFproj = regexprep(outputFolderCTFproj,'\.','p');
        mkdir(outputPath,outputFolderCTFproj);
        outputPathCTFproj   = [outputPath outputFolderCTFproj '/'];
    end
    inputPath = [inputParentPath DataFolderNames(kk).name '/'];
    loopcounter = 0;
    %% Loop over projections
    fprintf('Number of image processed and saved for data ''%s'':\n',DataFolderNames(kk).name)
    for nn =1:numel(dir([inputPath 'int_*.edf']))
        int = pmedfread(sprintf('%sint_%04u.edf',inputPath,nn));
        out = Reco(int,alpha,edp,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold);
        if evalTIElo
            edfwrite(sprintf('%sphase_%04u.edf',outputPathTIElo,nn),out.tieLO,'float32');
        end
        if evalTIEpnlo
            edfwrite(sprintf('%sphase_%04u.edf',outputPathTIEpnlo,nn),out.tiePNLO,'float32');
        end
        if evalCTF
            edfwrite(sprintf('%sphase_%04u.edf',outputPathCTF,nn),out.ctf,'float32');
        end
        if BinaryFilterThreshold(1) > 0
            edfwrite(sprintf('%sphase_%04u.edf',outputPathCTFproj,nn),out.ctfProjected,'float32');
        end
        loopcounter  = loopcounter + 1;
        totalcounter = totalcounter + 1;
        fprintf('%5u',nn)
        if mod(loopcounter,20) == 0
            fprintf('\n')
        end
    end
    fprintf('\n')
end
 telapsed = toc;
fprintf('A total of %u images processed and saved in %gs (%gs/image)\n',totalcounter,telapsed,telapsed/totalcounter)
