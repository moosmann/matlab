function DataProcID19(varargin)
% Data processing pipeline for experimental data taken at ID19@ESRF.
%
% Hot, dark, dead pixels filtering. Dark-and-flat-field correction. Interpolation of flat fields. Script
% assumes that the last 4 digits in filename are the angular position for
% flat field.
%
% ARGUMENTS
% ExpPath: string. Path to folder containing the subfolder 'data'
% DataSetPrefix: string. Default ''. Prefix specifying subset.
% OPTIONAL
% Experiment: string. Predefined set of parameters of a spedific Experiment
% DataSetsToProcess: integer vector. Default 0=ALL. Indices of data sets to process
%
% This script merges the old scripts FlatCor4cell (calling FlatCor) and
% RecoLoopXenopus4cell (calling RecoLoop) and includes additional processing
% steps: dark pixel filter, ring artifact filtering in sinogram,
% (tomographic reconsruction)
%
%Written by Julian Moosmann, 2014-01-17, last version: 2014-08-07

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%% PARSE INPUT ARGUMENT 'VARARGIN'
ParseVarargin(varargin)
%% DEFAULT ARGUMENTS AND PARAMETERS
if ~exist('Experiment','var')
    Experiment = '2014xeno17';
end
% Arguments set by individual parameters of specific experiment
switch lower(Experiment)
    case {'4cell','4-cell','xeno4cell','xeno4-cell','xenopus4cell','xenopus4-cell'}
        %ExpPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell';
        %ExpPath = '/export/scratch1/moosmann/ESRF_MI1079_ID19_July2011_inlineTomo';
        ExpPath = '/mnt/tomoraid-LSDF/users/moosmann/CWI_DATA/ESRF_MI1079_ID19_July2011_inlineTomo';
        DataSetPrefix = 'Xenopus';
        DataSetsToProcess = 1;
        projThresh = [1.024 0.95];%=[0.02 0.01];
        flatThresh = [0.02 0.001]; %[1.008,0.984];
        darkThresh = [0.02 0.001];%= [1.01, 0.985];
        EnergyDistancePixelsize = [20.0 0.945 0.75e-6];
    case {'zebra','zebrafish'}
    case {'2014xeno17'}
        ExpPath = '/mnt/LSDF/tomo/ESRF_ID19_2014-07-28_InHouse/';
        DataSetPrefix = '30jul04h03_stage17_e2v';
        DataSetsToProcess = 1;
        projThresh = [1.024 0.95];%=[0.02 0.01];
        flatThresh = [0.02 0.001]; %[1.008,0.984];
        darkThresh = [0.02 0.001];%= [1.01, 0.985];
        EnergyDistancePixelsize = [30.0 0.945 0.65e-6];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('DataSetPrefix','var')
    % Prefix of filename of projections. If empty
    DataSetPrefix = '';
end
if ~exist('filterSinoBeforePhase','var')
    % 0: none, 1: vertical (matlab), 2: horizontal (matlab), 3:  both
    filterSinoBeforePhase(1) = 0;
end
if ~exist('filterSinoAfterPhase','var')
    % 0: none, 1: vertical (matlab), 2: horizontal (matlab), 3:  both
    filterSinoAfterPhase(1) = 0;
end
if ~exist('maskingRatioThresh','var')
    % No masking for 0;
    maskingRatioThresh = 2;
end
% Renormalize intensity by its mean to yield I_fluc before phase retieval
if ~exist('renormInt','var')
    renormInt(1) = 1;
end
% Phase retrieval parameter
if ~exist('PhaseMethod','var')
    PhaseMethod = 'tie';'qp';'qpDual';'qp2';'tie';'ctf';
    % refractive index for H2O
    % Energy/eV  Delta,          Beta            Beta/Delta  -log10(Beta/Delta)
    % 20000      5.76620721E-07  3.46201373E-10  6.0040e-04  3.2216
    % 30000      2.56114134E-07  1.06431752E-10  4.1556e-04  3.3814
end
if ~exist('RegPar','var')
    % Standard PR: 2 seems superior to 2.5 for Xeno4cell data, 20 keV. RegPar = 1.5
    % still has prominent fringes and RegPar = 2.5 is contaminated by
    % large-scale variations
    %RegPar = 2;
    % Dual PR: RegPar = -log10(epsilon)
    RegPar = 2.5;2.5;3.22;
end
if ~exist('BinaryFilterThreshold','var')
    BinaryFilterThreshold = 0.1;
end
% Output formats: edf or tif
if ~exist('formatInt','var')
    formatInt = 'tif';
end
if ~exist('formatIntSino','var')
    formatIntSino = 'tif';
end
if ~exist('formatPhase','var')
    formatPhase = 'tif';
end
if ~exist('formatPhaseSino','var')
    formatPhaseSino = 'tif';
end
if ~exist('formatVolSlice','var')
    formatVolSlice = 'tif';
end
if ~exist('readInt','var')
    readInt(1) = 1;
end
% show figures, print output, etc
if ~exist('showFigures','var')
    showFigures = 0;
end
% Tomography
if ~exist('RotAxisPos','var')
   RotAxisPos = 1050; 
end
NumProjTomo = 1000;
FlatPrefix = 'ref';
%interpolation = 'nearest';'linear';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paths, folders, file names, etc
CheckTrailingSlash(ExpPath)
DataPath = [ExpPath 'data/'];
CheckTrailingSlash(DataPath);
IntPath = [ExpPath 'int/'];
CheckTrailingSlash(IntPath);
CheckAndMakePath(IntPath);
PhasePath = [ExpPath 'phase/'];
CheckTrailingSlash(PhasePath);
CheckAndMakePath(PhasePath);
VolPath = [ExpPath 'vol/'];
CheckTrailingSlash(VolPath);
CheckAndMakePath(VolPath);
DataFolderNames = FilenameCell([DataPath DataSetPrefix '*']);
% Data set indices to loop over
if isempty(DataSetsToProcess) || DataSetsToProcess(1)==0
    DataSetsToProcess = 1:numel(DataFolderNames);
end
switch filterSinoBeforePhase
    case {0,'none'}
        filterSinoBeforePhaseStr = '';
    case {1,'vertical'}
        filterSinoBeforePhaseStr = '_filtSinoVert';
    case {2,'horizontal'}
        filterSinoBeforePhaseStr = '_filtSinoHorz';
    case {3,'both'}
        filterSinoBeforePhaseStr = '_filtSinoBoth';
end
switch filterSinoAfterPhase
    case {0,'none'}
        filterSinoAfterPhaseStr = '';
    case {1,'vertical'}
        filterSinoAfterPhaseStr = '_filtSinoVert';
    case {2,'horizontal'}
        filterSinoAfterPhaseStr = '_filtSinoHorz';
    case {3,'both'}
        filterSinoAfterPhaseStr = '_filtSinoBoth';
end
hyphens(1:100) = '-';
fprintf('\n%s\nSTART DATA PROCESSING OF EXPERIMENT: %s\n',hyphens,Experiment);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dd = 1:numel(DataSetsToProcess)
    clear out;
    out.Experiment = Experiment;
    % Output path
    dataIndex = DataSetsToProcess(dd);
    out.DataSet = DataFolderNames{dataIndex};
    InputPath = [DataPath out.DataSet '/'];
    fprintf('\nDATA SET: %s\n',out.DataSet);
    fprintf('INPUT PATH: %s\n',InputPath);
    OutputPathInt     = sprintf('%s%s/int%s_%s',IntPath,out.DataSet,filterSinoBeforePhaseStr,formatInt);
    OutputPathIntSino = sprintf('%s%s/sino_int%s_%s',IntPath,out.DataSet,filterSinoBeforePhaseStr,formatIntSino);    
    % increment in loop, for testing purposes
    loopInc = 1;
    if ~readInt
        %% DARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(' Processing dark fields:\n');
        darkNames = FilenameCell([InputPath 'dark*']);
        if sum(strcmp('darkend0000.edf',darkNames))
            filename = [InputPath 'darkend0000.edf'];
            [dark, pixHot,pixDead,pixDark] = FilterPixel(pmedfread(filename)/15,darkThresh);
        elseif sum(strcmp('dark.edf',darkNames))
            filename = [InputPath 'dark.edf'];
            [dark, pixHot,pixDead,pixDark] = FilterPixel(pmedfread(filename),darkThresh);
        else
            dark = 0;
            fprintf('\n! NO DARK FIELD FOUND!\n');
        end
        out.dark.pixelsHot = pixHot;
        out.dark.pixelsDead = pixDead;
        out.dark.pixelsDark = pixDark;
        fprintf('  Elapsed time: %.2g s\n',toc);
        %% FLAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(' Processing flat fields:\n');
        % Check for HST files
        HSTNames  = FilenameCell([InputPath 'refHST*.edf']);
        NumRefPos  = numel(HSTNames);
        refMedian = zeros([NumRefPos size(dark)]);
        if  NumRefPos > 2
            % Loop over refHST files
            for nn = NumRefPos:-1:1
                % Read header and image of refHST file
                filename = [InputPath HSTNames{nn}];
                [h,im] = pmedf_read(filename);
                [im, pixHot,pixDead,pixDark] = FilterPixel(im,flatThresh);
                % Store additionall information: Flat-field postions, filtered
                % pixels, mean of refHST,ring current (if available)
                out.refHST.pixelsHot(nn) = pixHot;
                out.refHST.pixelsDead(nn) = pixDead;
                out.refHST.pixelsDark(nn) = pixDark;
                out.refHST.position(nn) = str2double(HSTNames{nn}(end-7:end-4));
                out.refHST.mean(nn) = mean(im(:));
                srcur = FindRingCurrentInEDFHeader(h);
                if ~isempty(srcur)
                    out.refHST.RingCurrent(nn) = srcur;
                end
                % Correct dark
                refMedian(nn,:,:) = im - dark;
            end
        else
            % Loop over different positions where flat fields were taken.
            RefPosStruct  = dir([InputPath  FlatPrefix '0000*.edf']); % Name struct of position where flat fields where recorded.
            NumRefPos = numel(RefPosStruct); % Number of flat-field positons.
            %RefPos         = struct('char',{},'num',{}); % Initialize a struct for the flat-field position.
            fcounter    = 0; % Counter of all recorded refs.
            pcounter    = 0; % Counter of flat-field positions.
            for nn = NumRefPos:-1:1
                pcounter     = pcounter + 1; % Total flat field counter.
                %RefPos(nn).char = RefPosStruct(nn).name(end-7:end-4); % Struct: String array of flat-field postions.
                %RefPos(nn).num  = str2double(RefPos(nn).char); % Struct: Number array of flat-field postions.
                out.refHST.position(nn) = str2double(RefPosStruct(nn).name(end-7:end-4));
                RefArray     = dir(sprintf('%s%s*_%04u.edf',InputPath,FlatPrefix,out.refHST.position(nn))); % Struct:
                NumRef       = numel(RefArray);
                RefCell      = cell(1,NumRef);
                RefHeader    = cell(1,NumRef);
                % Loop over flat fields taken at one position to improve statistics.
                for ll = NumRef:-1:1
                    fcounter = fcounter + 1;
                    [RefHeader{ll},RefCell{ll}] = pmedf_read([InputPath RefArray(ll).name]);
                end
                % Read ring current from edf header
                srcur = FindRingCurrentInEDFHeader(RefHeader{nn});
                if ~isempty(srcur)
                    out.refHST.ringCurrent(nn) = srcur;
                end
                refMedian(nn,:,:) = FilterPixel(median(cat(3,RefCell{:}),3)-dark,projThresh);
                % Mean of refHST
                %RefPos(nn).MeanOfRefMedian = mean(mean(refMedian(nn,:,:),3),2);
                % Save median filterd flat field.
                %edfwrite(sprintf('%smy_reffHST%s.edf',InputPath,RefPos(nn).char),refMedian(nn,:,:),'uint16');
                edfwrite(sprintf('%smy_reffHST%04u.edf',InputPath,out.refHST.position(nn)),squeeze(refMedian(nn,:,:)),'uint16');
                
            end
        end
        if maskingRatioThresh
            mask = MaskingAperture(dark,squeeze(min(refMedian,[],1)),maskingRatioThresh);
            out.maskingOfAperture.meanOfMask = mean(mask(:));
            out.maskingOfAperture.ratioThrehold = maskingRatioThresh;
        end
        fprintf('  Elapsed time: %.2g min\n',toc/60);
        %% PROJECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(' Read projections, filter pixels, flat-and-dark-field correction:\n');
        % Read file names
        ProjNames  = FilenameCell([InputPath DataSetPrefix '*.edf']);
        %NumProj  = numel(ProjNames);
        NumProj = NumProjTomo; %! only process projection within flat field range since the fast interpolation using interp1q does not extrapolate.
        stack = zeros([size(dark) NumProj]);
        %% Read, filter pixel, fd correct
        for nn = NumProj:-1*loopInc:1
            % Read
            PrintNum(nn);
            filename = [InputPath ProjNames{nn}];
            [h,im] = pmedf_read(filename);
            % Filter pixel
            [im, pixHot,pixDead,pixDark] = FilterPixel(im,projThresh);
            % Additional info: filtered pixels position, mean, ring current
            out.proj.position(nn) = nn;
            out.proj.pixelsHot(nn) = pixHot;
            out.proj.pixelsDead(nn) = pixDead;
            out.proj.pixelsDark(nn) = pixDark;
            out.proj.mean(nn) = mean(im(:));
            srcur = FindRingCurrentInEDFHeader(h);
            if ~isempty(srcur)
                out.proj.ringCurrent(nn) = srcur;
            end
            % Interpolated flat field correction
            %im = (im - dark)./squeeze(interp1(out.refHST.position,refMedian,nn,interpolation,'extrap'));
            im = (im - dark)./reshape( squeeze( interp1q(out.refHST.position,refMedian,nn) ), size(dark)  );
            % Aperture masking
            if maskingRatioThresh
                im = im.*mask;
                im = im + sum(im(:))/sum(mask(:))*(1-mask);
            end
            stack(:,:,nn) = im';
        end
    else
        fprintf(' Read intensities from disk:\n');
        % Read file names
        ProjNames  = FilenameCell([OutputPathInt '/*.' formatInt]);
        NumProj  = numel(ProjNames);
        filename = sprintf('%s/%s',OutputPathInt,ProjNames{1});
        stack = imread(filename,'tif');
        stack = zeros([size(stack) NumProj]);
        for nn = 1:loopInc:NumProj
            PrintNum(nn);
            filename = sprintf('%s/%s',OutputPathInt,ProjNames{nn});
            stack(:,:,nn) = imread(filename,'tif');
        end
    end
    fprintf('\n  Elapsed time: %.2g min\n',toc/60);
    %% Rotation axis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if RotAxisPos == 0
        dim1 = size(stack,1);
        y = round(dim1/2) + (-ceil(dim1/8):floor(dim1/8));
        RotAxisPos = ImageCorrelation( squeeze(stack(y,:,1)), fliplr(stack(y,:,end)) );
        RotAxisPos = RotAxisPos.HorizontalRotationAxisPosition;
    end
    %% Sino writing and optional filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n Filter sino before phase retrieval (%u) & write intensity sino:\n',filterSinoBeforePhase);
    CheckAndMakePath(OutputPathIntSino);
    for nn = 1:loopInc:size(stack,1)
        PrintNum(nn);
        % Filter sinogramm
        im = squeeze(stack(nn,:,:));
        if filterSinoBeforePhase
            im = FilterSino(im,2,filterSinoBeforePhase);
        end
        % Write sino slices
        filename = sprintf('%s/sino_int_%04u',OutputPathIntSino,nn);
        WriteImage(filename,im',formatIntSino);
    end
    fprintf('\n  Elapsed time: %.2g min\n',toc/60);
    
    %% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(' Write intensity maps & retrieve phase:\n');
    CheckAndMakePath(OutputPathInt);
    [phaseFilter,phaseAppendix] = PhaseFilter(PhaseMethod,size(stack),EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,'single');
    for nn = 1:loopInc:NumProj
        PrintNum(nn);
        im = squeeze(stack(:,:,nn));        
        % Write intensity maps
        if ~readInt
            filename = sprintf('%s/int_%04u',OutputPathInt,nn);
            WriteImage(filename,im,formatInt);
        end
        % Normalize intensity to one for phase retrieval
        if renormInt
            im = im/mean(im(:)) - 1;
            intStr = 'intFluc';
        else
            im = im - 1;
            intStr = 'int';
        end
        % Phase retrieval
        im = real(ifft2(phaseFilter.*fft2(im)));        
        stack(:,:,nn) = im;
    end
    fprintf('\n  Elapsed time: %.2g min\n',toc/60);    
    %% Filter sino, write phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n Filter sino after phase retrieval (%u) & write phase sino:\n',filterSinoAfterPhase);    
    OutputPathPhase = sprintf('%s%s/%s%s_phase%s_%s_%s',PhasePath,out.DataSet,intStr,filterSinoBeforePhaseStr,filterSinoAfterPhaseStr,phaseAppendix,formatPhase);
    OutputPathPhaseSino = sprintf('%s%s/sino_%s%s_phase%s_%s_%s',PhasePath,out.DataSet,intStr,filterSinoBeforePhaseStr,filterSinoAfterPhaseStr,phaseAppendix,formatPhaseSino);    
    OutputPathVol =              sprintf('%s%s/%s%s_phase%s_%s_%s',VolPath,out.DataSet,intStr,filterSinoBeforePhaseStr,filterSinoAfterPhaseStr,phaseAppendix,formatVolSlice);
    CheckAndMakePath(OutputPathPhaseSino);
    CheckAndMakePath(OutputPathVol);
    angles = 1*pi/(NumProjTomo) * ((1:NumProjTomo)-1);
    for nn = 1:loopInc:size(stack,1)
        PrintNum(nn);
        % Filter sinogramm
        im = squeeze(stack(nn,:,:));
        im = FilterSino(im,2,filterSinoAfterPhase);
        % Write sino
        filename = sprintf('%s/sino_phase_%04u',OutputPathPhaseSino,nn);
        WriteImage(filename,im',formatPhaseSino);
        % Reconstruct tomographic slice       
        sino = SubtractMean(RotAxisSymmetricCropping(im(:,1:NumProjTomo)',RotAxisPos));
        vol = astra_make_reco(sino,angles,'FBP_CUDA',1);
        % Write slice
        filename = sprintf('%s/slice_%04u',OutputPathVol,nn);        
        WriteImage(filename,vol,formatVolSlice);        
    end
    fprintf('\n  Elapsed time: %.2g min\n',toc/60);
    %% Write phase maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(' Write phase maps:\n');
    CheckAndMakePath(OutputPathPhase);    
    for nn = 1:loopInc:NumProj
        PrintNum(nn);
        filename = sprintf('%s/phase_%04u',OutputPathPhase,nn);
        WriteImage(filename,squeeze(stack(:,:,nn)),formatPhase);
    end
    fprintf('\n  Elapsed time: %.2g min\n',toc/60);    
    %% Processing finished.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.intensity.outputPath = OutputPathInt;
    out.intensity.sinoFiltering = filterSinoBeforePhase;
    out.phaseRetrieval.method = PhaseMethod;
    out.phaseRetrieval.regularizationParameter = RegPar;
    out.phaseRetrieval.outputPath = OutputPathPhase;
    out.phaseRetrieval.sinoFiltering = filterSinoAfterPhase;
    OutputPathMatlab = [ExpPath 'matlab/'];
    CheckTrailingSlash(OutputPathMatlab);
    CheckAndMakePath(OutputPathMatlab);    
    out.filename = sprintf('%s/%s_%s.mat',OutputPathMatlab,out.DataSet,datestr(now,'yyyy-mm-dd_HH-MM-SS'));
    out.timeElapsed = toc;
    save(out.filename,'out');
end
fprintf('\nTOTAL PROCESSING TIME: %g s = %.2g min = %.4g h\n%s\n',toc,toc/60,toc/3600,hyphens);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showFigures
    if ~readInt
        figure('Name','Mean of refHST VS postition')
        plot(out.refHST.position(:),out.refHST.mean(:))
        figure('Name','Mean of projection VS postition')
        plot(out.proj.position(:),out.proj.mean(:))
        figure('Name','Filtered hot pixel VS postition')
        plot(out.proj.position(:),out.proj.pixelsHot(:))
    end
end