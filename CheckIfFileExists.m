function CheckIfFileExists(FileNamePrefix,ParentPath)

%% Default arguments
if nargin < 1
    FileNamePrefix = '.txt';
end
if nargin < 2
    ParentPath = pwd;
end
%% Body
FolderStruct = dir(sprintf('%s/*',ParentPath));
FolderStruct(1:2) = [];
% Loop over folders
for n = 1:numel(FolderStruct)
    FileStruct = dir(sprintf('%s/*%s*',FolderStruct(n).name,FileNamePrefix));
    fprintf('%s: %u\n',FolderStruct(n).name,numel(FileStruct));
end