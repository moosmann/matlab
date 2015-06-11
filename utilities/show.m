function show(ImageNumber,FilePrefix,FilePostfix,WorkspaceVariableName)
% Reads tif or edf image in the current folder and assigns into a workspace
% variable. Reads also other image formats if defined in 'FilePostfix'
% argument. 
% !! If working directory contains thousands of files the function can be
% extremly slow since the (default) reading of all images contained in the
% current folder using Matlab function 'dir' takes ages!!
%
% ImageNumber: scalar>0, optional. if empty reads first tif image found the
% working directory by Matlab function 'dir'. default is 1.
% FilePrefix: string, optional. defines string pattern for images file names that
% should be read. default is [] so all tif or edf files are read.
% Affects 'ImageNumber'! 
% FilePostfix: string, optional. defines image file postfix: tif, edf, png,
% etc. default is [] then tif files are tried to read if fails edf files
% are tried. Affects 'ImageNumber'!   
% WorkspaceVariableName: string, optional. defines name of the workspace
% variable. default is [] then the first three letters of the image file
% name of the image are used.

%% Default arguments.
if nargin < 1
    ImageNumber = 1;
end
if nargin < 2
    FilePrefix = [];
end
if nargin < 3
    FilePostfix = [];
end
if nargin < 4
    WorkspaceVariableName = [];
end
%% Get image-file names in current folder
if isempty(FilePostfix)
    % Try tif files
    if ~isempty(FilePrefix) && ischar(FilePrefix)
        imStruct = dir([FilePrefix '*.tif']);
    else
        imStruct = dir('*.tif');
    end
    % Try edf files
    if isempty(imStruct)
        if ~isempty(FilePrefix) && ischar(FilePrefix)
            imStruct = dir([FilePrefix '*.edf']);
        else
            imStruct = dir('*.edf');
        end
    end
    % Assign FilePostfix
    FilePostfix = imStruct(1).name(end-2:end);
else
    if ~isempty(FilePrefix) && ischar(FilePrefix)
        imStruct = dir([FilePrefix '*.' FilePostfix]);
    else
        imStruct = dir(['*.' FilePostfix]);
    end
end
imName = imStruct(ImageNumber).name;
%% Read image
switch FilePostfix
    case 'tif'
        im = double(imread(imName));
    case 'edf'
        im = pmedfread(imName)';
end
%% Show image
figure('Name',imName)
imshow(im,[],'InitialMagnification','fit')
colorbar
%% Assign image to workspace
if isempty(WorkspaceVariableName)
   WorkspaceVariableName = [imName(1:3) num2str(ImageNumber,'%04u')];
end
if ischar(WorkspaceVariableName)
    assignin('base',WorkspaceVariableName,im);
end
%% Print info about image
fprintf('Read %s image ''%s'' of size [%4u x %4u] and assigned it to workspace variable ''%s''.\n', ...
    FilePostfix,imName,size(im),WorkspaceVariableName)