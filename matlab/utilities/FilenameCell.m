function CellOfFilenames = FilenameCell(DirectoryName)
% Output is a cell containing the filenames found in the directory
% 'DirectoryName'. Default directoy is present working directory. Output
% cell does not include '.' and '..'.
%
% Written by Julian Moosmann

if nargin < 1
    DirectoryName = './*';
end

CellOfFilenames = dir(DirectoryName);

% for aLoop = 1:numel(CellOfFilenames)
%     if length(strfind(CellOfFilenames{aLoop},'ref')) > 1
%         CellOfFilenames{aLoop} = {};
%     end
% end

CellOfFilenames = {CellOfFilenames(:).name};

if ~isempty(CellOfFilenames)
    while CellOfFilenames{1}(1) == '.'
        CellOfFilenames(1) = [];
    end
end