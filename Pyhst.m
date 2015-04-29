function Pyhst(PathToParFiles,ParFilePrefix,DataSetsToRecon)
% Function reading 'par' files in given folder and calling PyHST from
% within MATLAB. If DataSetsToRecon is empty (default), the function gives
% a prompt and waits for input. If DataSetsToRecon is 0, all DataSets found
% are reconstructed. If DataSetsToRecon is a 1x?-vector the corresponding
% data sets are reconstructed.

%% Default arguments.
if nargin < 1
    % Path to the folder containing the par file.s
   PathToParFiles = pwd;
   %'/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/vol/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';
end
if nargin < 2
    % Part of the prefix of the names of the par files:
    % sprintf('*%s*.vol.par',ParFilePrefix)
    ParFilePrefix = '';
end
if nargin < 3
    % [] (empty): user is prompted to input the data sets that shall be
    % reconstructed, scalar: 0=all data sets will be reconstructed,
    % 1x?-row-vector of the data sets to be reconstructed
    DataSetsToRecon = [];
end
%% Body.
% Check ending of path string.
if PathToParFiles(end) ~= '/'
    PathToParFiles = [PathToParFiles '/'];
end
% Read par files and print names to screen.
ParFileStruct = dir([PathToParFiles '*' ParFilePrefix '*.vol.par']);
NumParFiles = numel(ParFileStruct);
fprintf('\nPAR FILES FOUND WITH PREFIX ''%s'' IN DIRECTORY\n%s\n',ParFilePrefix,PathToParFiles)
for nn = 1:NumParFiles
    fprintf(' %2u: %s\n',nn,ParFileStruct(nn).name)
end
if isempty(DataSetsToRecon)
    % Call input function.
    DataSetNum = input(sprintf('TYPE IN AS A ROW VECTOR THE NUMBER OF THE TOMOGRAMs TO RECONSTRUCT [1..%u]: ',NumParFiles));
else
    if DataSetsToRecon == 0;
        DataSetNum = 1:NumParFiles;
    else
        DataSetNum = DataSetsToRecon;
    end
end
% Loop over tomo set to be reconstructed.
if (min(DataSetNum(:)) > 0) && (max(DataSetNum(:)) <= NumParFiles)
    % Print par file names that will be reconstructed
    fprintf('DATA SETs THAT WILL BE RECONSTRUCTED:\n')
    for nn = 1:numel(DataSetNum)
        fprintf(' %2u: %s\n',DataSetNum(nn),ParFileStruct(DataSetNum(nn)).name)
    end
    % Call PyHST within MATLAB
    fprintf('START LOOPING PyHST OVER PAR FILES.\n')
    for nn = 1:numel(DataSetNum)
        unix(sprintf('for f in %s%s; do echo $f; /opt/pyhst/PyHST $f; done',PathToParFiles,ParFileStruct(DataSetNum(nn)).name));%/opt/pyhst/PyHST $f
    end
else
    fprintf('NO TOMOGRAMS WILL BE RECONSTRUCTED.\n')
end
