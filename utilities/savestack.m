function savestack( vol, save_path, im_type)
% Save image stack as image sequence.
%
% vol : 3d image stack
% save_path : absolute or relative path
% im_type : string. default: tif
% 

if nargin < 3
    im_type = 'tif';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = toc;

CheckAndMakePath( save_path);
CheckTrailingSlash( save_path );

fprintf( '\nSaving volume [%u %u %u] to %s', size( vol), save_path )
switch im_type
    case 'tif'
        parfor nn = 1:size( vol, 3)
            filename = sprintf('%sreco_%04u.tif', save_path, nn);
            im = vol(:,:,nn);            
            write32bitTIF( filename, im);
        end
end
fprintf( '\nDone in %.1f s (%.2f min).\n', toc -t, (toc - t)/60 )