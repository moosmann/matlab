function etool(ImageNumber,FilenamePattern,WorkspaceVariableName)
% Read edf file and display image.
%
% Written by Julian Moosmann, last version 2013-1024

%% Default arguments
if nargin < 1
    ImageNumber = 0;
end
if nargin < 2
    FilenamePattern = '';
end
if nargin < 3
    WorkspaceVariableName = '';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read edf file names in current folder
if isempty(FilenamePattern)
    imStruct = FilenameCell('*.edf');
else
    imStruct = FilenameCell([FilenamePattern '.edf']);
end
NumIm = numel(FilenameCell);
%% Set image number
if ImageNumber == 0
    ImageNumber = ceil(NumIm/2);
end
%% Read image
imName = imStruct{ImageNumber};
im = pmedfread(imName)';
%% Show image
itool(im,imName)
%% Assign image to workspace
if isempty(WorkspaceVariableName)
    if ~isempty(regexp(imName(1),'\d', 'once'))
        WorkspaceVariableName = lower(['im' imName(1:3) num2str(ImageNumber,'%04u')]);          
    else
        WorkspaceVariableName = lower([imName(1:3) num2str(ImageNumber,'%04u')]);          
    end   
end
if ischar(WorkspaceVariableName)
    assignin('base',WorkspaceVariableName,im);
end
%% Print info
fprintf('Read ''edf'' image ''%s'' of %u. Size: %4u x %4u. Variable ''%s'' is assigned to base workspace.\n',imName,NumIm,size(im),WorkspaceVariableName)
