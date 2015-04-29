function Pipeline32ID201307()
% Script which (interactively) loops over the scans of the experiment In
% vivo imaging, GUP at beamline 32ID-C at APS, ANL, Chicago, Illinois, USA,
% 2013-07-28 to 2013-08-02.
%
%Written Julian Moosmann, last modified 2013-09-19

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/';
%% Write matlab output into logfile
LogFileFolder = '/mnt/tomoraid-LSDF/users/moosmann/matlabLog/';
CheckTrailingSlash(LogFileFolder);
[~, hostname] = unix('echo $HOSTNAME');
hostname = hostname(1:end-1);
LogFile = sprintf('%s%s_%s',LogFileFolder,hostname,datestr(now,'yyyy-mm-dd_HH-MM-SS'));
diary(LogFile); % Start logging of Matlab output using the function diary
try
    fprintf('\nSCRIPT: %s',mfilename);
    fprintf('\nSTART: %s',datestr(now,'yyyy.mm.dd  HH:MM:SS'));
    fprintf('\nHOSTNAME: %s',hostname);
    fprintf('\nLOG FILE: %s',LogFile);
    %% Interactive query
    DataSets = ChooseDataSet([ParentPath 'data']);
    TomoToProcess  = input(sprintf('\nTOMOGRAMS TO RECONSTRUCT ([INDEX AS IT APPEARS IN FILENAME], DEFAULT=ALL): '));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Loop over different stages/embryos
    for nn = 1:numel(DataSets)
        try
            fprintf('\nSTART RECONSTRUCTION OF DATA SET NUMBER %2u: %s\n',nn,DataSets{nn});
            %% Start loop over time-lapse tomograms
            DataSetName = DataSets{nn};        
            DataProc32IDstack201307tomo('DataSet',DataSetName,'TomoToProcess',TomoToProcess,'PixelRegion',{},'saveIntMaps','edf','savePhaseMaps','edf');
        catch err
            report = getReport(err);
            fprintf('\n\nERROR/EXCEPTION OCCURED WHILE LOOPING OVER DATA SETS:\n\n%s\n\n',report);
        end
    end
    matlabpool close
    diary off
catch err
    report = getReport(err);
    fprintf('\n\nERROR/EXCEPTION OCCURED:\n\n%s\n\n',report);
    diary off
    matlabpool close
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataSetsForReco = ChooseDataSet(DataPath)
% Interactive command line prompt to choose the data sets to be
% reconstructed
%% Default arguments
if nargin < 1
    DataPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/data';
end
DataSetsForReco = {};
%% Read PARENT data set names into struct
CheckTrailingSlash(DataPath);
ParDataSets = FilenameCell([DataPath '*']);
while ParDataSets{1}(1) == '.'
    ParDataSets(1) = [];
end
NumParDataSets  = numel(ParDataSets);
% Print found parent data sets 
fprintf('\nDATA SETS FOUND IN ''%s'':\n',DataPath)
for nn = 1:NumParDataSets
    fprintf(' %2u: %s\n',nn,ParDataSets{nn})
end
% Prompt for PARENT data sets to be reconstructed
DataSetInd = input(sprintf('CHOOSE PARENT FOLDER OF DATA SETS ([1-%u],DEFAUT=ALL): ',NumParDataSets));
if ~isempty(DataSetInd)
    ParDataSets = ParDataSets(DataSetInd);
end
NumParDataSets  = numel(ParDataSets);
%% Prompt if all data sets found under parent folder shoud be reonstructed
RecoAll = input(sprintf('RECONSTRUCT ALL DATA SETS FOUND IN PARENT FOLDERS (0=NO=DEFAULT,1=YES)? '));
if isempty(RecoAll)
    RecoAll = 0;
end
%% Read data set names into output cell
for nn = 1:NumParDataSets
    DataSets = FilenameCell([DataPath ParDataSets{nn} '/*']);
    while DataSets{1}(1) == '.'
        DataSets(1) = [];
    end
    NumDataSets = numel(DataSets);
    fprintf('\nDATA SETS FOUND IN ''%s'':\n',ParDataSets{nn})
    for mm = 1:NumDataSets
        fprintf(' %2u: %s\n',mm,DataSets{mm});
        DataSets{mm} = [ParDataSets{nn} '/' DataSets{mm}];
    end
    % Prompt for data sets to be reconstructed
    if RecoAll
        DataSetInd = [];
    else
        DataSetInd = input(sprintf('CHOOSE DATA SETS ([1-%u],DEFAUT=ALL): ',NumDataSets));
    end
    if ~isempty(DataSetInd)
        DataSets = DataSets(DataSetInd);
    end
    DataSetsForReco = cat(2,DataSetsForReco,DataSets);
end
