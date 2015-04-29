function eshow(ImageNumber,FilenamePattern,WorkspaceVariableName)
% Read edf file, display image, and assign it to workspace.

%% Default arguments
if nargin < 1
    ImageNumber = 1;
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
%% Read image
imName = imStruct{ImageNumber};
im     = pmedfread(imName)';
%% Show image
ishow(im,imName)
%% Assign image to workspace
if isempty(WorkspaceVariableName)
   WorkspaceVariableName = lower([imName(1:3) num2str(ImageNumber,'%04u')]);
end
if ischar(WorkspaceVariableName)
    assignin('base',WorkspaceVariableName,im);
end
%% Print info
fprintf('Read ''edf'' image ''%s'' of %u. Size: %4u x %4u. Variable ''%s'' is assigned to base workspace.\n',imName,NumIm,size(im),WorkspaceVariableName)
