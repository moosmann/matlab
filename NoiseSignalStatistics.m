function NoiseSignalStatistics(ParentPath,FlatPrefix,DataSetsToProcess,ROIfactor,HotPixThreshold,NumProjAv,doShow)
% Compute Poisson noise of the flat fields and statistics of the (flat
% field corrected) intensity data.
%
% This file is the combinatin of the two functions 'NoiseStatistics' and
% 'PhaseContrastSignal'. The former is computational extensive, therefore
% the latter was written after the former to be seperately executable.
% However due to bad data (overexposed flat fields in two sets) and a bug in the former
% function both functions are now combined into one.

%% Default arguments.
if nargin < 1
    ParentPath = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo/';
end
if nargin < 2
    FlatPrefix = 'ref';
end
if nargin < 3
    DataSetsToProcess = 0;
end
if nargin < 4
    %first row: horizontal ROI, secod row: vertical ROI, corresponds to
    %Matlab's matrix notation;
    ROIfactor = 4;
end
if nargin < 5
    HotPixThreshold = 0.01;
end
if nargin < 6
    NumProjAv = 40;
end
if nargin < 7
    doShow = 0;
end

fprintf('\nNOISE AND SIGNAL STATISTICS OF FLAT FIELDS AND SAMPLE DATA\n')
%% Name string of path, folder and files.
% Default setting of 'ParentPath'.
if ParentPath == 0
    ParentPath = pwd;
end
% Chechk ending of string.
if ParentPath(end) == '/'
    ParentPath(end) = [];
end
% 0000000001111111111222222222233333333333444444
% 1234567890123456789012345678901234567890123456
% opt_xeno_stage11_tomo_22keV_10x_125mm_0p20sec_
%% Open txt for writing.
% Text file for verbose output of flat field analysis.
txtFileNameVerbose = [ParentPath '/NoiseSignalStatisticsVerboseOutput.txt'];
fidVer = fopen(txtFileNameVerbose,'wt');
fprintf(fidVer,'POISSON STATISTICS OF FLAT FIELDS\n');
fprintf('Parent path of data sets: %s\n',ParentPath)
fprintf(fidVer,'Parent path of data sets: %s\n',ParentPath);
fprintf('Percentage of hot-filtered pixels in flat fields: %g\n',HotPixThreshold)
fprintf(fidVer,'Percentage of hot-filtered pixels in flat fields: %g\n',HotPixThreshold);
fprintf('ROI factor for flat fields: x = %g (the first and the last xth are cut off)\n',ROIfactor)
fprintf(fidVer,'ROI factor for flat fields: x = %g (the first and the last xth are cut off)\n',ROIfactor);
fprintf('ROI factor for sample: x = %g (a xth from the image is centrically cropped)\n',ROIfactor)
% Text file in table format.
txtFileNameTable = [ParentPath '/NoiseSignalStatisticsTableFormat.txt'];
fidTab = fopen(txtFileNameTable,'wt');
fprintf(fidTab,'stage E z dx dt <sam> samStd samStd*<sam> samStd/<sam> samMin samMax samMax-samMin <flatCounts> <flatVar> <flatConvFac> <photons> <%%PoissonNoise>\n');
fprintf(fidTab,'- keV mm micron ms 1 1 1 1 1 1 1 1 1 1 1 1\n');
% Print to standard output.
fprintf('Text file for verbose output: %s\n',txtFileNameVerbose)
fprintf('Text file for table format: %s\n',txtFileNameTable)
%% Define intput paths and folder names.
% Set path to the folder where the data sets are contained in.
DataPath = [ParentPath '/data/'];
IntPath = [ParentPath '/int/'];
% Default prefix of the names of the folders the data is in.
% Get data folder names.
DataFolderStruct      = dir(DataPath);
IntFolderStruct       = dir(IntPath);
DataFolderStruct(1:2) = [];
IntFolderStruct(1:2)  = [];
if DataSetsToProcess == 0
    DataSetsToProcess = 1:numel(DataFolderStruct);
end
NumOfSets           = numel(DataSetsToProcess); 
out(NumOfSets).mean = 0;
out(NumOfSets).std  = 0;
DistanceOffset = 13; %mm
% 1/CropFac*Dim will be cropped center symmetrically.
CropFac = 3;
%% Loop over data sets.
for jj = NumOfSets:-1:1
    fprintf('\nPROCESSING DATA SET: %s\n',DataFolderStruct(DataSetsToProcess(jj)).name);
    %% SIGNALS ANALYSIS
    IntFolderName = IntFolderStruct(jj).name;
    SubIntPath       = [IntPath IntFolderName '/'];
    %% Read parameters from directory names.
    % 0000000001111111111222222222233333333333444444
    % 1234567890123456789012345678901234567890123456
    % opt_xeno_stage11_tomo_22keV_10x_125mm_0p20sec_
    out(jj).Stage     =  str2double(IntFolderName(15:16));
    out(jj).Energy    =  str2double(IntFolderName(23:24));%in keV
    out(jj).Distance  = (DistanceOffset + str2double(IntFolderName(33:35)));%in mm
    out(jj).Pixelsize = 15/str2double(IntFolderName(29:30));%in microns
    out(jj).ExpoTime  = str2double(IntFolderName(41:42))*10;%in milliseconds
    %fprintf('%2u. %-50s: ',jj,IntFolderName);
    fprintf(' stage:         %2u\n energy:        %2u keV\n pixel size:   %4.02g micron\n distance:      %3u mm\n exposure time: %3u ms',...
        out(jj).Stage,out(jj).Energy,out(jj).Pixelsize,out(jj).Distance,out(jj).ExpoTime)
   
    %% Calculate mean values and standard deviations.
    fs       = dir([SubIntPath 'int*.edf']);
    NumFiles = length(fs)-6;
    for ff = NumProjAv:-1:0
        % Read image and get dimension.
        im = pmedfread([SubIntPath fs(floor(NumFiles/NumProjAv*ff)+1).name])';
        [dimx dimy] = size(im);
        % Crop central region.
        xx = ceil(dimx/2*(CropFac-1)/CropFac):floor(dimx/2*(CropFac+1)/CropFac);
        yy = ceil(dimy/2*(CropFac-1)/CropFac):floor(dimy/2*(CropFac+1)/CropFac);
        im = im(xx,yy);
        imMean(ff+1) = mean(im(:));
        imStd(ff+1)  = std(im(:));
        imMax(ff+1)  = max(im(:));
        imMin(ff+1)  = min(im(:));
    end
    % Average of mean values and standard deviations.
    out(jj).mean     = mean(imMean(:));
    out(jj).std      = mean(imStd(:));
    out(jj).stdmean  = out(jj).std*out(jj).mean;
    out(jj).std2mean = out(jj).std/out(jj).mean;
    out(jj).max      = mean(imMax(:));
    out(jj).min      = mean(imMin(:));
    fprintf(fidTab,'%u %u %u %.2f %u %g %g %g %g %g %g %g ',...
        out(jj).Stage,out(jj).Energy,out(jj).Distance,out(jj).Pixelsize,out(jj).ExpoTime,out(jj).mean,out(jj).std,out(jj).stdmean,out(jj).std2mean,out(jj).min,out(jj).max,out(jj).max-out(jj).min);
    fprintf('\n SAMPLE: MEAN:%5.3g',out(jj).mean)
    fprintf('\n         STD:%6.3g',out(jj).std) 
    fprintf(' STD*MEAN:%7.3g',out(jj).stdmean)
    fprintf(' STD/MEAN:%6.3g',out(jj).std2mean)
    fprintf('\n         MIN:%5.3g',out(jj).min) 
    fprintf(' MAX:%5.3g',out(jj).max)
    fprintf(' MAX-MIN:%5.3g\n',out(jj).max-out(jj).min)
    % Print more info.
    if doShow
        ishow(im)
        fprintf('        MEANs: %s\n',mat2str(imMean,3))
        fprintf('        STDs: %s\n',mat2str(imStd,3))
    end
    %% FLAT FIELD ANALYSIS
    % Loop over different positions where flat fields were taken.
    SubDataPath = [DataPath DataFolderStruct(DataSetsToProcess(jj)).name '/'];
    % Region of interest. !!!!! edf images are not transposed after !!!!!
    if ROIfactor > 0
        im = pmedfread([SubDataPath 'ref0000_0000.edf']);
        [dim1 dim2] = size(im);
        ROIhorizontal = floor(dim1/ROIfactor):ceil(dim1*(1-1/ROIfactor));
        ROIvertical   = floor(dim2/ROIfactor):ceil(dim2*(1-1/ROIfactor));
    end
    fprintf(fidVer,'\nData set: %s\n',DataFolderStruct(DataSetsToProcess(jj)).name);
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
    fprintf(' FLAT: 2D means of 1D means along 3rd dim: %s\n',mat2str(RefMeanMean,Precision))
    fprintf('       2D means of 1D vars along 3rd dim:  %s\n',mat2str(RefVarMean,Precision))
    fprintf('       conversion factors:                 %s\n',mat2str(PropFac,Precision))
    fprintf('       numbers of photons:                 %s\n',mat2str(NumPhotons,Precision+1))
    fprintf('       Poisson noise in %%:                 %s\n',mat2str(PoissonNoise,Precision))
    % Write to text file.
    fprintf(fidVer,'2D means of 1D means along 3rd dim: %s\n',mat2str(RefMeanMean,Precision));
    fprintf(fidVer,'2D means of 1D vars along 3rd dim:  %s\n',mat2str(RefVarMean,Precision));
    fprintf(fidVer,'conversion factors:                 %s\n',mat2str(PropFac,Precision+1));
    fprintf(fidVer,'numbers of photons:                 %s\n',mat2str(NumPhotons,Precision+1));
    fprintf(fidVer,'Poisson noise in %%:                %s\n',mat2str(PoissonNoise,Precision));
    % Mean over flat field positions.
    RefMeanMeanMean(jj) = mean(RefMeanMean(:));
    RefVarMeanMean(jj)  = mean(RefVarMean(:));
    PropFacMean(jj)     = mean(PropFac(:));
    NumPhotonsMean(jj)  = mean(NumPhotons(:));
    PoissonNoiseMean(jj)= mean(PoissonNoise(:));
    % Print to screen.
    fprintf('       <detector counts>:     %g\n',RefMeanMeanMean(jj))
    fprintf('       <variance>:            %g\n',RefVarMeanMean(jj))
    fprintf('       <conversion factor>:   %g\n',PropFacMean(jj))
    fprintf('       <photons at detector>: %g\n',NumPhotonsMean(jj))
    fprintf('       <Poisson noise in %%>:  %g\n',PoissonNoiseMean(jj))
    % Write to text file.
    fprintf(fidVer,'<detector counts>:     %g\n',RefMeanMeanMean(jj));
    fprintf(fidVer,'<variance>:            %g\n',RefVarMeanMean(jj));
    fprintf(fidVer,'<conversion factor:    %g\n',PropFacMean(jj));
    fprintf(fidVer,'<photons at detector>: %g\n',NumPhotonsMean(jj));
    fprintf(fidVer,'<Poisson noise in %%>:  %g\n',PoissonNoiseMean(jj));
    fprintf(fidTab,' %.1f %.1f %.6f %.1f %.6f\n',...
        RefMeanMeanMean(jj),RefVarMeanMean(jj),PropFacMean(jj),NumPhotonsMean(jj),PoissonNoiseMean(jj));
    % Clear.
    clear RefMeanMean RefVarMean PropFac NumPhotons PoissonNoise
end
fclose(fidVer);
fclose(fidTab);
