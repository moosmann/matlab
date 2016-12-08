function CellOfFilenames = FilenameCell(DirectoryName)
% Output is a cell containing the filenames found in the directory
% 'DirectoryName'. Default directoy is present working directory. Output
% cell does not include '.' and '..'.
%
% Written by Julian Moosmann, last modified 2013-09-13

if nargin < 1
    DirectoryName = './*';
end

CellOfFilenames = dir(DirectoryName);
CellOfFilenames = {CellOfFilenames(:).name};

if ~isempty(CellOfFilenames)
    while CellOfFilenames{1}(1) == '.'
        CellOfFilenames(1) = [];
    end
end