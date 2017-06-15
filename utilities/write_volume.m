function write_volume( tag, vol, output_type, reco_path, raw_bin, reco_bin, counter_offset, verbose )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
    counter_offset = 0;
end
if nargin < 8
    verbose = 1;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tag == 0
    return
else
    t = toc;
    switch output_type
        
        %% Single precision: 32-bit float tiff
        case 'float'
            
            if reco_bin <= 1
                % Unbinned
                PrintVerbose(verbose, '\n Write floats:')
                save_path = sprintf( '%sfloat_rawBin%u/', reco_path, raw_bin);
            else
                % Binned
                PrintVerbose(verbose, '\n Write floats binned:')
                save_path = sprintf( '%sfloat_rawBin%u_recoBin%u/', reco_path, raw_bin, reco_bin);
            end
            
            CheckAndMakePath( save_path )
            
            parfor mm = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_path, mm + counter_offset);
                write32bitTIFfromSingle( filename, vol( :, :, mm) )
            end
            
            %% 16-bit tiff
        case 'uint16'
            
            t = toc;
            if reco_bin <= 1
                %% Unbinned
                PrintVerbose(verbose, '\n Write uint16:')
                save_path = sprintf( '%s/uint16_rawBin%u/', reco_path, raw_bin);
            else
                %% Binned
                PrintVerbose(verbose, '\n Write uint16 binned:')
                save_path = sprintf( '%s/uint16_rawBin%u_recoBin%u/', reco_path, raw_bin, reco_bin);
            end
            CheckAndMakePath( save_path )
            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_path, nn + counter_offset);
                imwrite( uint16( (2^16 - 1) * vol( :, :, nn) ), filename );
            end
            
            %% 8-bit tiff
        case 'uint8'
            
            if reco_bin <= 1
                % Unbinned
                PrintVerbose(verbose, '\n Write uint8:')
                save_path = sprintf( '%s/uint8_rawBin%u/', reco_path, raw_bin);
                
            else
                % Binned
                PrintVerbose(verbose, '\n Write uint8 binned:')
                save_path = sprintf( '%s/uint8_rawBin%u_recoBin%u/', reco_path, raw_bin, reco_bin);
            end
            
            CheckAndMakePath( save_path )
            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_path, nn + counter_offset);
                imwrite( uint8( (2^8 - 1) * vol( :, :, nn) ), filename );
            end
            
    end
    PrintVerbose(verbose, ' done in %.2f min.', (toc - t) / 60)
end
