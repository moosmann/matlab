function FlatCorThis(doWrite,InputPath,DataPrefix,FlatPrefix,OutputPath,RegionOfInterest,HotPixFiltThreshold,HotPixFiltPrintInfo,doReadRingCurrent)
% Dark-and-flat-field correction of tomography data, including
% interpolation of flat fields. Script assumes that the last 4 digits in
% filename are the angular position for flat field.
% 
% Written by J. Moosman
% First version: 21/09/2011
% Last modification: 21/09/2011 (if updated)
%
% doWrite: signed integer (default: 0). For doWrite > 0 the processed
% images are saved, for < 0 the images are only processe without saving
% (for testing). For doWrite < 0 and > 0, every x-th projection is
% processed where x is the modulus of doWriten. Images are saved in:
% ParentPath/int/Scan%/, see next argument.
% ParentPath: string or 0 (default). Parent path of the data set, where
% following data structure is assumed: 
% ParentPath/data/Scan%/%.edf. If 0 the current working directory
% (pwd) is used. "/" at the end of the string is mandatory.
% DataPrefix: sting or 0 (default). Part of the prefix of data file name. Doesn't need to be
% complete prefix, only has to be unambigous. If 0 first four letter of the
% name of the folder (Scane%) where the data is contained in were taken.
% FlatPrefix: string. Part of the prefix of flat-field file name. Doesn't need to be
% complete, only has to be unambigous. Default is 'ref'.
% RegionOfInterest: cell of 2x2 matrices (for each data set one hast to
% define the ROI), 2x2 matrix (ROI is apllied to each data set), or 0
% (default, whole image is processed). Matrix of start and end pixels of
% the region of interest. If 0 the whole image is processed. 
% HotPixFiltThreshold: scalar (0.01, default). If < 1 then interpreted as
% the ratio of all pixel which are to be filtered. If > 1 it's the
% threshold used to determine which pixel should be filtered (for details
% see function 'FilterHotPixel'.
% HotPixFiltPrintInfo: 0 or 1 (default). Print output of the function
% 'FilterHotPixel'.
% doReadRingCurrent: 1 or 0 (default). Reads the ring current contained in
% the header of the flat-field images if available, alternatively calculates
% the mean of each flat-field image and then show a plot of it. Meant for
% detecting jumps in the photon flux which complicate the procedure of the
% flat-field correction due to interpolation (not yet implemented)

%% Default arguments.
if nargin < 1
    doWrite = 0;
end
if nargin < 2
    % String: trailing seperator ('/') not needed. E.g.
    InputPath = pwd;
end
if nargin < 3
    DataPrefix = 'xeno';
end
if nargin < 4
    FlatPrefix = 'ref';
end
if nargin < 5
    OutputPath = [pwd '/int/'];
end
if nargin < 6
    %first row: horizontal ROI, secod row: vertical ROI, corresponds to
    %Matlab's matrix notation;
    RegionOfInterest = 0;
end
if nargin < 7
    HotPixFiltThreshold = 0.03;
end
if nargin < 8
    HotPixFiltPrintInfo = 1;
end
if nargin < 9
    doReadRingCurrent = 0;
end

%% Check input and output paths.
if InputPath(end) ~= '/'
    InputPath = [InputPath '/'];
end
outputPath = OutputPath;
if outputPath(end) ~= '/'
    outputPath = [outputPath '/'];
end
if ~exist(outputPath,'dir')
    mkdir(outputPath)
end

%% Region of interest. !!!!! edf images are not transposed after !!!!!
if ~isscalar(RegionOfInterest) && ~iscell(RegionOfInterest)
    ROIhorizontal = RegionOfInterest(1,1):RegionOfInterest(1,2);
    ROIvertical = RegionOfInterest(2,1):RegionOfInterest(2,2);
end
fprintf('\nApply flat-/dark-field correction including hot-pixel filtering\n');
fprintf('PART OF PREFIX OF DATA FILE NAMES: %s\n',DataPrefix);
fprintf('INPUT PATH:  %s\n',InputPath);
fprintf('OUTPUT PATH: %s\n',outputPath);
%% Read dark field if existing.
fprintf('PROCESS DARK FIELD:\n');
if HotPixFiltPrintInfo
    fprintf('  ');
end
darkStruct = dir([InputPath 'darkend0000.edf']);
if numel(darkStruct) == 1
    im = pmedfread([InputPath 'darkend0000.edf'])/15;
    if ~isscalar(RegionOfInterest)
        dark = im(ROIhorizontal,ROIvertical);
    else
        dark = im;
    end
    dark = FilterHotPixel(dark,HotPixFiltThreshold,HotPixFiltPrintInfo);
else
    dark = 0;
end
%% Check for HST files.
fprintf('PROCESS FLAT FIELDs:\n');
hstPosArray  = dir([InputPath 'refHST*.edf']);
NumOfRefPos  = numel(hstPosArray);
if  NumOfRefPos > 2
    fprintf('  refHST images found.\n');
    if doReadRingCurrent
        fprintf('  [FlatPosition     < int > RingCurrent] = \n');
    end
    for kk = NumOfRefPos:-1:1
        [RefHeader{kk} im] = pmedf_read([InputPath hstPosArray(kk).name]);
        if ~isscalar(RegionOfInterest)
            RefMedian{kk} = im(ROIhorizontal,ROIvertical);
        else
            RefMedian{kk} = im;
        end
        if HotPixFiltPrintInfo
            fprintf('  ');
        end
        RefMedian{kk} = FilterHotPixel(RefMedian{kk}-dark,HotPixFiltThreshold,HotPixFiltPrintInfo);
        RefPos(kk).char = hstPosArray(kk).name(end-7:end-4); % Struct: String array of flat-field postions.
        RefPos(kk).num  = str2double(RefPos(kk).char); % Struct: Number array of flat-field postions.
        RefMedianMean(kk) = mean(RefMedian{kk}(:));
        if doReadRingCurrent
            % Read ring current. !!!! Modification needed in line 74 in
            % function "pmedf_findPos". "omega_pos = sscanf(motors_pos, '%g','C')" -> "omega_pos = sscanf(motors_pos, '%g')"; !!!!
            counters.mne = pmedf_findInHeader(RefHeader{kk}, 'counter_mne', 'string');
            if ~isempty(counters.mne)
                counters.pos = pmedf_findInHeader(RefHeader{kk}, 'counter_pos', 'string');
                srcur.pos    = pmedf_findMnePos(counters.mne, 'srcur');
                srcur.value  = pmedf_findPos(counters.pos, srcur.pos);
                fprintf('  [%12u%12.5g%12.5g]\n',kk,RefMedianMean(kk),srcur.value);
            else
                fprintf('  [%12u%12.5g]\n',kk,RefMedianMean(kk));
            end
        end
    end
    if doReadRingCurrent
        figure('Name','mean of refHST image VS projection index')
        plot(squeeze(cat(3,RefPos.num)),RefMedianMean,'-+')
    end
else
    %% Loop over different positions where flat fields were taken.
    RefPosArray  = dir([InputPath  FlatPrefix '0000*.edf']); % Name struct of position where flat fields where recorded.
    NumOfRefPos = numel(RefPosArray); % Number of flat-field positons.
    RefPos         = struct('char',{},'num',{}); % Initialize a struct for the flat-field position.
    %RefMedian   = cell(1,NumOfRefPos); % Initialize cell for the median filtered flat fields.
    tread       = 0; % Read time.
    tmedian     = 0; % Time for median filtering.
    fcounter    = 0; % Counter of all recorded refs.
    pcounter    = 0; % Counter of flat-field positions.
    %fprintf(1,'(ProjIndex RefIndex <int>) = \n');
    for kk = NumOfRefPos:-1:1
        pcounter     = pcounter + 1; % Total flat field counter.
        RefPos(kk).char = RefPosArray(kk).name(end-7:end-4); % Struct: String array of flat-field postions.
        RefPos(kk).num  = str2double(RefPos(kk).char); % Struct: Number array of flat-field postions.
        RefArray     = dir([InputPath FlatPrefix '*_' RefPos(kk).char '.edf']); % Struct:
        NumRef       = numel(RefArray);
        RefCell      = cell(1,NumRef);
        RefHeader    = cell(1,NumRef);
        % Loop over flat fields taken at one position to improve statistics.
        for ll = NumRef:-1:1
            fcounter = fcounter + 1;
            tic
            [RefHeader{ll} im] = pmedf_read([InputPath RefArray(ll).name]);
            if ~isscalar(RegionOfInterest)
                RefCell{ll} = im(ROIhorizontal,ROIvertical);
            else
                RefCell{ll} = im;
            end
            tread = tread + toc;
        end
        if doReadRingCurrent
            % Read ring current. !!!! Modification needed in line 74 in
            % function "pmedf_findPos". "omega_pos = sscanf(motors_pos, '%g','C')" -> "omega_pos = sscanf(motors_pos, '%g')"; !!!!
            counters.mne = pmedf_findInHeader(RefHeader{ll}, 'counter_mne', 'string');
            counters.pos = pmedf_findInHeader(RefHeader{ll}, 'counter_pos', 'string');
            srcur.pos    = pmedf_findMnePos(counters.mne, 'srcur');
            srcur.value  = pmedf_findPos(counters.pos, srcur.pos);
            fprintf(1,'  [%12u%12.5g%12.5g]\n',kk,mean(RefCell{ll}(:)),srcur.value);
        end
        % Compute median of flat fields at given flat-field position.
        tic
        RefMedian{kk} = FilterHotPixel(median(cat(3,RefCell{:}),3)-dark,HotPixFiltThreshold,HotPixFiltPrintInfo);
        %RefMedian{kk} = sum(cat(3,RefCell{:}),3)/NumRef;
        %domain(RefMedian{kk},'median filtered flat field')
        RefMedianMean(kk) = mean(RefMedian{kk}(:));
        % Save median filterd flat field.
        tmedian = tmedian + toc;
        if doWrite > 0
            % !! If detection of refHST files fails, they will be
            % overwritten !!
            edfwrite(sprintf('%s%sHST%s.edf',InputPath,FlatPrefix,RefPos(kk).char),RefMedian{kk},'uint16');
        end
    end
    fprintf(1,'Read %2u flats in %fs. Computed median filtered flat fields at %u angular positions each of %2u images in %fs\n',fcounter,tread,pcounter,fcounter/pcounter,tmedian);
    clear RefCell
    if doReadRingCurrent
        figure('Name','mean of median filtered reference image VS projection index')
        plot(squeeze(cat(3,RefPos.num)),RefMedianMean,'-+')
    end
end
%% Interpolate flat fields and make flat-and-dark field correction.
DataPosArray  = dir([InputPath DataPrefix '*.edf']); % Name struct of position where data was taken.
NumOfDatPos  = numel(DataPosArray); % Number of angular positions.
%DatPos       = struct('char',{},'num',{}); % Initialize a struct for the data positions.
icounter     = 0; % Counter for sample images
tic;
%% Loop over projections.
fprintf(1,'PROCESSING projection number: \n');
HotPixFiltPrintInfoForInts = 0;
% Decrease number of projection to process for testing (doWrite=0).
if doWrite == 1
    ProjNumIncrement = 1;
elseif doWrite > 1
    ProjNumIncrement = doWrite;
elseif doWrite == 0
    ProjNumIncrement = 400;
elseif doWrite < 0
    ProjNumIncrement = abs(doWrite);
end
for nn = 1:ProjNumIncrement:NumOfDatPos
    fprintf(1,' %4u',nn);
    icounter = icounter + 1;
    if mod(icounter,20)==0
        fprintf(1,'\n');
    end
    im = pmedfread([InputPath DataPosArray(nn).name]);
    if ~isscalar(RegionOfInterest)
        int = im(ROIhorizontal,ROIvertical);
    else
        int = im;
    end
    int = (int - dark)./squeeze(interp1(shiftdim(cat(3,RefPos.num),1),shiftdim(cat(3,RefMedian{:}),2),nn,'linear','extrap'));
    int = FilterHotPixel(int,HotPixFiltThreshold,HotPixFiltPrintInfoForInts);
    if doWrite > 0
        edfwrite(sprintf('%s/int_%04u.edf',outputPath,nn),int,'float32');
    end
end;
telpased = toc;
fprintf(1,'\nInterpolation of reference images and flat-/dark-field correction done at %u points in %gs (%gs/image)\n\n',icounter,telpased,telpased/icounter);
%% Assign variables in user's current workspace.
if doWrite <= 0
    %assignin('base','RefPosArray',RefPosArray);
    assignin('base','RefPos',RefPos);
    assignin('base','RefMedian',RefMedian);
    assignin('base','RefMedianMean',RefMedianMean);
    assignin('base','dark',dark);
    assignin('base','DataPosArray',DataPosArray);
    assignin('base','int',int);
end
