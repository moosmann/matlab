function write_volume_tif( vol, save_path, im_class, im_prefix)
% Save image stack as image sequence.
%
% vol : 3D image stack
% save_path : absolute path
% im_class : image class, 'float', 'uint16', or 'unit8', if empty input
%   class of 'vol' is used. if integer class is used, 'vol' is assumed to
%   be in [0,1].
% im_prefix: filename prefix to be used
%
% Writte by J. Moosmann

if nargin < 3
    im_class = 'float';
end
if nargin < 4
    im_prefix = 'im';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = toc;

CheckAndMakePath( save_path);
CheckTrailingSlash( save_path );

if isempty(im_class)
    im_class = class(vol);
end
convert = ~strcmp(class(vol), im_class);

fprintf( '\nSaving volume [%u %u %u] to \n %s', size( vol), save_path )
switch im_class
    case 'float'
        parfor nn = 1:size( vol, 3)
            filename = sprintf('%s%s_%06u.tif', save_path, im_prefix, nn);            
            write32bitTIFfromSingle( filename, vol(:,:,nn) )
        end
        
    case 'uint16'
        parfor nn = 1:size( vol, 3)
            filename = sprintf('%s%s_%06u.tif', save_path, im_prefix, nn);
            if convert
                imwrite( uint16( (2^16 - 1) * vol(:,:,nn) ), filename );
            else
                imwrite( vol(:,:,nn), filename );
            end
        end
        
    case 'uint8'
        parfor nn = 1:size( vol, 3)
            filename = sprintf('%s%s_%06u.tif', save_path, im_prefix, nn);
            if convert
                imwrite( uint8( (2^8 - 1) * vol(:,:,nn) ), filename );
            else
                imwrite( vol(:,:,nn), filename );
            end
        end
        
    otherwise
        error( 'Class ''%s'' not supported.', output_type )
end
fprintf( '\n done in %.1f s (%.2f min).\n', toc -t, (toc - t)/60 )