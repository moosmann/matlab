function save_path = write_volume( tag, vol, output_type, reco_path, raw_bin, phase_bin, reco_bin, counter_offset, verbose, suffix )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 8
    counter_offset = 0;
end
if nargin < 9
    verbose = 1;
end
if nargin < 10
    suffix = '';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tag == 0
    return
else
    t = toc;
    
    save_path = sprintf( '%s%s_rawBin%u', reco_path, output_type, raw_bin);
    if phase_bin > 1
        save_path = sprintf( '%s_phaseBin%u', save_path, phase_bin);
    end
    PrintVerbose(verbose, '\n Write %s', output_type)
    if reco_bin > 1
        PrintVerbose(verbose, ' binned')
        save_path = sprintf( '%s_recoBin%u', save_path, reco_bin);
    end
    PrintVerbose(verbose, ':')
    save_path = [save_path suffix];
    CheckTrailingSlash( save_path )
    CheckAndMakePath( save_path )
        
    switch output_type
        
        case 'float'            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_path, nn + counter_offset);
                write32bitTIFfromSingle( filename, vol( :, :, nn) )
            end
            
        case 'uint16'            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_path, nn + counter_offset);
                imwrite( uint16( (2^16 - 1) * vol( :, :, nn) ), filename );
            end
            
        case 'uint8'            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_path, nn + counter_offset);
                imwrite( uint8( (2^8 - 1) * vol( :, :, nn) ), filename );
            end
        otherwise
           error( 'Format ''%s'' not supported.', output_type )
    end
    PrintVerbose(verbose, ' done in %.2f min.', (toc - t) / 60)
end
