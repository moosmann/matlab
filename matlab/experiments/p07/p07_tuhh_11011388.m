%TUHH B5

%% Extract scan names to be stitched
fprintf( '\n\nIDENTIFY SCANS TO BE STITCHED:' )
proc_path = sprintf('/asap3/petra3/gpfs/p07/2021/data/11011388/processed/');
d = dir( proc_path );
s = {};
sn = 1;
for n = 1:numel(d)
    name1 = d(n).name;
    if strcmp( name1(end), 'a')
        name2 = d(n+1).name;
        if strcmp( name2(end), 'b')
            name0 = name1(1:end-2);            
            s{sn} = name0;            
            fprintf( '\n %4u %s %s %s', sn, name1, name2, name0)
            sn = sn + 1;
            
        end
    end
end

%% Stitch scans
fprintf( '\n\nSTITCH SCANS' )
for n = 135%1:numel(s) 
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