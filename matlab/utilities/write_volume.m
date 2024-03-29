function save_path = write_volume( tag, vol, output_type, write, raw_bin, phase_bin, reco_bin, counter_offset, verbose, suffix)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 8
    counter_offset = 0;
end
if nargin < 9
    verbose = 0;
end
if nargin < 10
    suffix = '';
end

if ~write.is_phase
    reco_path = assign_from_struct( write, 'reco_path', 0 );
else
    reco_path = assign_from_struct( write, 'reco_phase_path', 0 );
end
deleteFiles = assign_from_struct( write, 'deleteFiles', 0 );
beamtimeID_regexp  = assign_from_struct( write, 'beamtimeID', '' );
scan_name = assign_from_struct( write, 'scan_name', '' );

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tag == 0
    return
else
    save_path = sprintf( '%s%s_rawBin%u', reco_path, output_type, raw_bin);
    if phase_bin > 1
        save_path = sprintf( '%s_phaseBin%u', save_path, phase_bin);
    end
    PrintVerbose(verbose, '\nWrite %s', output_type)
    if reco_bin > 1
        PrintVerbose(verbose, ' binned')
        save_path = sprintf( '%s_recoBin%u', save_path, reco_bin);
    end
    PrintVerbose(verbose, ':')
    save_path = [save_path suffix];
    CheckTrailingSlash( save_path )
    CheckAndMakePath( save_path, deleteFiles, beamtimeID_regexp )
    
    if strcmp( output_type(1:4), 'uint' )
        tlow = write.tlow;
        thigh = write.thigh;
        vol = -tlow + vol;
        vol = 1 / (thigh - tlow) * vol;
    end
    switch output_type
        
        case 'float'
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%s_%06u.tif', save_path, scan_name, nn - 1 + counter_offset);
                write32bitTIFfromSingle( filename, vol( :, :, nn) )
            end
            
        case 'float_adapthisteq'
            imah = @(im) (adapthisteq(normat(im)));
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%s_%06u.tif', save_path, scan_name, nn - 1 + counter_offset);
                im = vol(:,:,nn);
                im = imah(im);
                write32bitTIFfromSingle(filename,im)
            end
            
        case 'uint16'
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%s_%06u.tif', save_path, scan_name, nn - 1  + counter_offset);
                imwrite( uint16( (2^16 - 1) * vol( :, :, nn) ), filename );
            end
            
        case 'uint8'            
            parfor nn = 1:size( vol, 3)
                filename = sprintf( '%sreco_%s_%06u.tif', save_path, scan_name, nn - 1 + counter_offset);
                imwrite( uint8( (2^8 - 1) * vol( :, :, nn) ), filename );
            end
        otherwise
           error( 'Format ''%s'' not supported.', output_type )
    end
end
