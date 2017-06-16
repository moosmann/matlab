function tshow(ImageNumber,FilenamePattern,WorkspaceVariableName, workspace)
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
if nargin < 4
    workspace = 'caller';'base';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read edf file names in current folder
if isempty(FilenamePattern)
    imStruct = FilenameCell('*.tif');
else
    imStruct = FilenameCell([FilenamePattern '.tif']);
end
NumIm = numel(imStruct);
%% Set image number
if ~ImageNumber
    ImageNumber = ceil(NumIm/2);
end
if ImageNumber < 0
    ImageNumber = NumIm + ImageNumber + 1;
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
    assignin(workspace, WorkspaceVariableName, im);
end
%% Print info
fprintf('Read ''%s'' of %u tiff images foundin current directory.\nImage shape: %4u %4u.\nVariable ''%s'' is assigned to %s workspace.\n',imName,NumIm,size(im),WorkspaceVariableName, workspace)
