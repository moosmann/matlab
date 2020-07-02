function NoiseStatistics(ParentPath,FlatPrefix,DataSetsToProcess,ROIfactor,HotPixThreshold)
% Compute Poisson noise from flat fields.

%% Default arguments.
if nargin < 1
    ParentPath = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo/';
end
if nargin < 2
    FlatPrefix = 'ref';
end
if nargin < 3
    DataSetsToProcess = 29;
end
if nargin < 4
    %first row: horizontal ROI, secod row: vertical ROI, corresponds to
    %Matlab's matrix notation;
    ROIfactor = 4;
end
if nargin < 5
    HotPixThreshold = 0.00;
end
%% Name string of path, folder and files.
% Default setting of 'ParentPath'.
if ParentPath == 0
    ParentPath = pwd;
end
% Chechk ending of string.
if ParentPath(end) == '/'
    ParentPath(end) = [];
end
% Text file name string.
txtFileName = [ParentPath '/PoissonStatisticsOfFlatFields2.txt'];
fprintf('Text file: %s\n',txtFileName)
fid = fopen(txtFileName,'wt');
% Print.
fprintf('\nPOISSON STATISTICS OF FLAT FIELDS\n')
fprintf(fid,'POISSON STATISTICS OF FLAT FIELDS\n');
fprintf('Parent path of data sets: %s\n',ParentPath)
fprintf(fid,'Parent path of data sets: %s\n',ParentPath);
fprintf('Percent of hot-filtered pixels: %g\n',HotPixThreshold)
fprintf(fid,'Percent of hot-filtered pixels: %g\n',HotPixThreshold);
fprintf('ROI factor: x = %g (the first and the last xth are cut off)\n',ROIfactor)
fprintf(fid,'ROI factor: x = %g (the first and the last xth are cut off)\n',ROIfactor);
% Set path to the folder where the data sets are contained in.
DataPath = [ParentPath '/data/'];
% Default prefix of the names of the folders the data is in.
% Get data folder names.
DataFolderNames      = dir(DataPath);
DataFolderNames(1:2) = [];
if DataSetsToProcess == 0
    DataSetsToProcess = 1:numel(DataFolderNames);
end
%% Loop over data sets.
for jj = numel(DataSetsToProcess):-1:1
    % Loop over different positions where flat fields were taken.
    SubDataPath = [DataPath DataFolderNames(DataSetsToProcess(jj)).name '/'];
    % Region of interest. !!!!! edf images are not transposed after !!!!!
    if ROIfactor > 0
        im = pmedfread([SubDataPath 'ref0000_0000.edf']);
        [dim1 dim2] = size(im);
        ROIhorizontal = floor(dim1/ROIfactor):ceil(dim1*(1-1/ROIfactor));
        ROIvertical   = floor(dim2/ROIfactor):ceil(dim2*(1-1/ROIfactor));
    end
    fprintf('\nPROCESSING DATA SET: %s\n',DataFolderNames(DataSetsToProcess(jj)).name);
    fprintf(fid,'\nData set: %s\n',DataFolderNames(DataSetsToProcess(jj)).name);
    RefPosArray = dir([SubDataPath  FlatPrefix '0000*.edf']); % Name struct of position where flat fields where recorded.
    NumOfRefPos = numel(RefPosArray); % Number of flat-field positons.
    RefPos      = struct('char',{},'num',{}); % Initialize a struct for the flat-field position.
    % Read dark field.
    darkStruct = dir([SubDataPath 'darkend0000.edf']);
    if numel(darkStruct) == 1
        im = pmedfread([SubDataPath 'darkend0000.edf'])/15;
        if ROIfactor > 0
            dark = im(ROIhorizontal,ROIvertical);
        else
            dark = im;
        end
        if HotPixThreshold > 0
            dark = FilterHotPixel(dark,HotPixThreshold,0);
        end
    else
        dark = 0;
    end
    for kk = NumOfRefPos:-1:1
        RefPos(kk).char = RefPosArray(kk).name(end-7:end-4); % Struct: String array of flat-field postions.
        RefPos(kk).num  = str2double(RefPos(kk).char); % Struct: Number array of flat-field postions.
        RefArray     = dir([SubDataPath FlatPrefix '*_' RefPos(kk).char '.edf']); % Struct:
        NumRef       = numel(RefArray);
        RefCell      = cell(1,NumRef);
        RefHeader    = cell(1,NumRef);
        % Loop over flat fields taken at one position to improve statistics.
        for ll = NumRef:-1:1
            [RefHeader{ll} im] = pmedf_read([SubDataPath RefArray(ll).name]);
            if HotPixThreshold > 0
                im                 = FilterHotPixel(im,HotPixThreshold,0);
            end
            if ROIfactor > 0
                RefCell{ll} = im(ROIhorizontal,ROIvertical);
            else
                RefCell{ll} = im;
            end
        end
        % 1-D mean of different flat fields at position 'kk'.
        RefMean{kk}     = mean(cat(3,RefCell{:})-repmat(dark,[1 1 NumRef]),3);
        RefVar{kk}      = var( cat(3,RefCell{:}),1,3);
        % 2-D mean of 1-D mean of flat fields at position 'kk'.
        RefMeanMean(kk)  = mean(RefMean{kk}(:));
        RefVarMean(kk)   = mean(RefVar{kk}(:));
        PropFac(kk)      = RefVarMean(kk)/RefMeanMean(kk);
        NumPhotons(kk)   = RefMeanMean(kk)./PropFac(kk);
        PoissonNoise(kk) = 100*(NumPhotons(kk)).^(-1/2);
    end
    % Print to screen.
    Precision = 5;
    fprintf('2D mean of 1D mean along 3rd dim: %s\n',mat2str(RefMeanMean,Precision))
    fprintf('2D mean of 1D var along 3rd dim:  %s\n',mat2str(RefVarMean,Precision))
    fprintf('Multiplication factor: %s\n',mat2str(PropFac,Precision))
    fprintf('Number of Photons: %s\n',mat2str(NumPhotons))
    fprintf('Poisson noise in %%: %s\n',mat2str(PoissonNoise,Precision))
    % Write to text file.
    fprintf(fid,'2D mean of 1D mean along 3rd dim: %s\n',mat2str(RefMeanMean,Precision));
    fprintf(fid,'2D mean of 1D var along 3rd dim:  %s\n',mat2str(RefVarMean,Precision));
    fprintf(fid,'Multiplication factor: %s\n',mat2str(PropFac,Precision+1));
    fprintf(fid,'Number of Photons: %s\n',mat2str(NumPhotons));
    fprintf(fid,'Poisson noise in %%: %s\n',mat2str(PoissonNoise,Precision));
    % Mean over flat field positions.
    RefMeanMeanMean(jj) = mean(RefMeanMean(:));
    RefVarMeanMean(jj)  = mean(RefVarMean(:));
    PropFacMean(jj)     = mean(PropFac(:));
    NumPhotonsMean(jj)  = mean(NumPhotons(:));
    PoissonNoiseMean(jj)= mean(PoissonNoise(:));
    % Print to screen.
    fprintf('Mean detector counts: %g\n',RefMeanMeanMean(jj))
    fprintf('Mean variance: %g\n',RefVarMeanMean(jj))
    fprintf('Mean proportionality factor: %g\n',PropFacMean(jj))
    fprintf('Mean real photons: %g\n',NumPhotonsMean(jj))
    fprintf('Mean Poisson noise: %g\n',PoissonNoiseMean(jj))
    % Write to text file.
    fprintf(fid,'Mean detector counts: %g\n',RefMeanMeanMean(jj));
    fprintf(fid,'Mean variance: %g\n',RefVarMeanMean(jj));
    fprintf(fid,'Mean proportionality factor: %g\n',PropFacMean(jj));
    fprintf(fid,'Mean real photons: %g\n',NumPhotonsMean(jj));
    fprintf(fid,'Mean Poisson noise: %g\n',PoissonNoiseMean(jj));
    % Clear.
    clear RefMeanMean RefVarMean PropFac NumPhotons PoissonNoise
end
fclose(fid);
