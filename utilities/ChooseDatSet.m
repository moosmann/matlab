function FolderStruct = ChooseDatSet(PathToDataSets,DataSetPrefix)

%% Default arguments
if nargin < 1
    PathToDataSets = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/data/';
end
if nargin < 2
    DataSetPrefix = '';
end
%% Body
% Check path string
if PathToDataSets(end) ~= '/'
    PathToDataSets = [PathToDataSets '/'];
end
% Read data set names into struct.
if isempty(DataSetPrefix)
    FolderStruct = dir(PathToDataSets);
    FolderStruct(1:2) = [];
else
    FolderStruct = dir([PathToDataSets '*' DataSetPrefix '*']);
end
NumDataSets  = numel(FolderStruct);
% Print names of data sets found in 'PathToDataSets'
fprintf('\nDATA SETS FOUND IN THE DIRECTORY %s\n',PathToDataSets)
for nn = 1:NumDataSets
    fprintf(' %2u: %s\n',nn,FolderStruct(nn).name)
end
% Prompt for data set(s) that should be reconstructed
DataSetNum = input(sprintf('ENTER AS ROW VECTOR THE NUMBER OF DATA SETs TO RECONSTRUCT [1-%u]: ',NumDataSets));
if ~isempty(DataSetNum)
    FolderStruct = FolderStruct(DataSetNum);
end
%fprintf('DATA SET WHICH WILL BE RECONSTRUCTED: %s\n',mat2str(DataSetNum))