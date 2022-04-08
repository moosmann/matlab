function WriteVol(vol, path, filename_prefix, slicing_dimension)
% Loop over input volume and save slice.



%% Defaults
if nargin < 2
    path = pwd;
end
if nargin < 3
    filename_prefix = 'slice';
end
if nargin < 4
    slicing_dimension = 3;
end

%% Main
% make dir
% path = sprintf('%s_%04ux%04ux%04u',path,size(vol));
% if nargin > 1
%     CheckAndMakePath(path);
% end

% loop
for nn = 1:size(vol,slicing_dimension)
    switch slicing_dimension
        case 3
            filename = sprintf('%s/%s%07u',path,filename_prefix,nn);
            WriteImage(filename,squeeze(vol(:,:,nn)),'tif')
        case 2
            filename = sprintf('%s/%s%07u',path,filename_prefix,nn);
            WriteImage(filename,squeeze(vol(:,nn,:)),'tif')
        case 1
            filename = sprintf('%s/%s%07u',path,filename_prefix,nn);            
            WriteImage(filename,squeeze(vol(nn,:,:)),'tif')
    end
end