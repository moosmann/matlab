%TUHH B5

for n = 14:20
    %scan_path = '/asap3/petra3/gpfs/p07/2021/data/11011388/processed/tuhh_003_000_0p55mm';
    scan_path = sprintf('/asap3/petra3/gpfs/p07/2021/data/11011388/processed/tuhh_003_%03u_0p55mm', n );
    scan_subfolder = 'reco_phase/tie_regPar2p00';
    reco_subfolder = 'float_rawBin2';
    crop = 1;
    save_stitched_volume = 1;
    stitched_volume_path = '';
    
    %[s, vol] = stitch_volumes( scan_path, scan_subfolder, reco_subfolder, crop, save_stitched_volume, stitched_volume_path );
    stitch_volumes( scan_path, scan_subfolder, reco_subfolder, crop, save_stitched_volume, stitched_volume_path );
    
    
    
end