function [name_cell, full_path_cell] = FilenameCell(full_path)
% Output is a cell containing the filenames found in the directory
% 'full_path'. Default directoy is present working directory. Output
% cell does not include '.' and '..'.
%
% Written by Julian Moosmann

if nargin < 1
    full_path = './*';
end

d = dir(full_path);

% for aLoop = 1:numel(name_cell)
%     if length(strfind(name_cell{aLoop},'ref')) > 1
%         name_cell{aLoop} = {};
%     end
% end

name_cell = {d(:).name};
full_path_cell = cellfun( @(a,b) [a filesep b], {d(:).folder}, {d(:).name}, 'UniformOutput', 0 );

if ~isempty(name_cell)
    while name_cell{1}(1) == '.'
        name_cell(1) = [];
        full_path_cell(1) = [];
    end
end

