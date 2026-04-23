p = '/asap3/petra3/gpfs/p07/2026/data/11022041/processed/';
scans = {
'004_gst_no1_c0_'
'005_gst_no2_c0_'
'006_gst_No1_c1A_'
'007_gst_No2_c1d_'
'008_gst_No1_c1d_'
'009_gst_No2_c1d_'
'010_gst_No3_c27A_'
'011_gst_No1_c3A_'
'012_gst_No2_c3A_'
'013_gst_No1_c3D_'
'014_gst_No2_c3D_'
'015_gst_No4_c28a_'
'016_gst_No1_c5a_'
'017_gst_No2_c5a_'
'018_gst_No1_c5d_'
'019_gst_No2_c5d_'
'020_gst_No3_c27d_'
'021_gst_No4_c27d_'
'022_gst_No5_c27A_'
'023_gst_No6_c27A_'
'024_gst_No5_c27d_'
'025_gst_No7_c27a_'
'026_gst_No6_c27d_'
'027_gst_No7_c27d_'
};
scan_subfolder = 'reco';
reco_subfolder = 'float_rawBin2';
stitched_volume_path = [];
m1 = [0 1 1 1 1 1 1 1 1 1 1];
m2 = [];%[1 1 1 1 1 1 1 1 1 1 1];
for n = 24:numel(scans)
    scan_path = [p scans{n}];
    if n < 19
        scan_mask = m1;
        fprintf('\n %2u: %s %u',n,scan_path,scan_mask(1))
    else
        scan_mask = m2;
        fprintf('\n %2u: %s no mask',n,scan_path)
    end
    s = stitch_volumes(scan_path,scan_subfolder,reco_subfolder,stitched_volume_path,scan_mask);
end
fprintf('\n FINISHED')