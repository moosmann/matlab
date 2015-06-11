function WriteVol(vol,FolderName,FilenamePrefix,LoopDirection)
% Loop over input volume and save slice.



%% Defaults
if nargin < 2
    FolderName = 'vol';
end
if nargin < 3
    FilenamePrefix = 'slice';
end
if nargin < 3
    LoopDirection = 3;
end

%% Main
% make dir
FolderName = sprintf('%s_%04ux%04ux%04u',FolderName,size(vol));
if nargin > 1
    CheckAndMakePath(FolderName);
end

% loop
for nn = 1:size(vol,LoopDirection)
    switch LoopDirection
        case 3
            filename = sprintf('%s/%s_%04u',FolderName,FilenamePrefix,nn);
            WriteImage(filename,squeeze(vol(:,:,nn)),'tif')
        case 2
            filename = sprintf('%s/%s_%04u',FolderName,FilenamePrefix,nn);
            WriteImage(filename,squeeze(vol(:,nn,:)),'tif')
        case 1
            filename = sprintf('%s/%s_%04u',FolderName,FilenamePrefix,nn);            
            WriteImage(filename,squeeze(vol(nn,:,:)),'tif')
    end
end