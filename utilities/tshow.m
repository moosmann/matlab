function tshow(ImageNumber,FilenamePattern,WorkspaceVariableName)
% Read tif file and display image.
%
% Written by Julian Moosmann, last version 2013-1024

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    ImageNumber = false;
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
    imStruct = FilenameCell('*.tif');
else
    imStruct = FilenameCell([FilenamePattern '.tif']);
end
NumIm = numel(FilenameCell);
%% Set image number
if ~ImageNumber
    ImageNumber = ceil(NumIm/2);
end
%% Read image
imName = imStruct{ImageNumber};
im = double(imread(imName));
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
fprintf('Read ''tif'' image ''%s'' of %u. Size: %4u x %4u. Variable ''%s'' is assigned to base workspace.\n',imName,NumIm,size(im),WorkspaceVariableName)





