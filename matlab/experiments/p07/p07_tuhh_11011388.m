%TUHH B5

%% Extract scan names to be stitched
% parent path
proc_path = sprintf('/asap3/petra3/gpfs/p07/2021/data/11011388/processed/');
d = dir( proc_path ); % struct with all folders/files in parent path
s = {}; % cell to store folder names to be stitched
sn = 1; % counter
fprintf( '\n\nIDENTIFY SCANS TO BE STITCHED:' )
fprintf( '\n %4s  % s %s  %s', '#', '1st_height', '2nd_height', 'sitched_name')
% loop over all folders found and identify scan names ending with 'a' and 'b'
for n = 1:numel(d)
    name1 = d(n).name;
    % identfify scan names ending with 'a'
    if strcmp( name1(end), 'a')
        name2 = d(n+1).name;
        % check if the following scan ends with 'b'
        if strcmp( name2(end), 'b')
            name0 = name1(1:end-2);
            % store common scan name part in cell
            s{sn} = name0;            
            fprintf( '\n %4u  %s  %s  %s', sn, name1, name2, name0)
            sn = sn + 1;
        end
    end
end
fprintf( '\n' )

%% Stitch scans
fprintf( '\nSTITCH SCANS' )
for n = 1:numel(s) 
    name = s{n};
    scan_path = sprintf( '%s%s', proc_path, name);
    scan_subfolder = 'reco_phase/tie_regPar2p00';
    reco_subfolder = 'float_rawBin2';
    crop = 1;
    save_stitched_volume = 1;    
    stitched_volume_path = '';
    
    full_path_a = sprintf( '%s_a/%s/%s', scan_path, scan_subfolder, reco_subfolder );
    full_path_b = sprintf( '%s_b/%s/%s', scan_path, scan_subfolder, reco_subfolder );
    
    fprintf( '\n\n %4u %s %u %u', n, scan_path, exist( full_path_a, 'dir' ), exist( full_path_b, 'dir' ) )
    
    %[s, vol] = stitch_volumes( scan_path, scan_subfolder, reco_subfolder, crop, save_stitched_volume, stitched_volume_path );
    stitch_volumes( scan_path, scan_subfolder, reco_subfolder, crop, save_stitched_volume, stitched_volume_path );
end
fprintf( '\n' )