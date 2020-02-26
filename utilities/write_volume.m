function save_path = write_volume( tag, vol, output_type, par, reco_bin, counter_offset, suffix)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 8
    counter_offset = 0;
end
if nargin < 9
    par.verbose = 0;
end
if nargin < 10
    suffix = '';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~par.phase_retrieval
    reco_path = par.write_reco_path;
else
    reco_path = par.write_reco_phase_path;
end
scan_name = par.scan_name; 


if tag == 0
    return
else
    save_path = sprintf( '%s%s_rawBin%u', reco_path, output_type, par.raw_bin);
    if par.phase_retrieval_post_binning_factor > 1
        save_path = sprintf( '%s_phaseBin%u', save_path, par.phase_bin);
    end
    PrintVerbose(par.verbose, '\n Write %s', output_type)
    if reco_bin > 1
        PrintVerbose(par.verbose, ' binned')
        save_path = sprintf( '%s_recoBin%u', save_path, reco_bin);
    end
    PrintVerbose(par.verbose, ':')
    save_path = [save_path suffix];
    CheckTrailingSlash( save_path )
    CheckAndMakePath( save_path, par.write_deleteFiles, par.write_beamtimeID_regexp )
        
    switch output_type
        
        case 'float'            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%s_%06u.tif', save_path, scan_name, nn + counter_offset);
                write32bitTIFfromSingle( filename, vol( :, :, nn) )
            end
            
        case 'uint16'            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%s_%06u.tif', save_path, scan_name, nn + counter_offset);
                imwrite( uint16( (2^16 - 1) * vol( :, :, nn) ), filename );
            end
            
        case 'uint8'            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%s_%06u.tif', save_path, scan_name, nn + counter_offset);
                imwrite( uint8( (2^8 - 1) * vol( :, :, nn) ), filename );
            end
        otherwise
           error( 'Format ''%s'' not supported.', output_type )
    end
end
