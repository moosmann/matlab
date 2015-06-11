dataPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/data/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/';
savePath = '/mnt/tomoraid-LSDF/users/moosmann/NatureProtocols/preprocessing/images/';
proj51 = double(imread([dataPath 'proj_00051.tif']));
proj51_hpf = FilterHotPixel(proj51,0.02,1);
write32bitTIF([savePath 'Proj51.tif'],proj51);
write32bitTIF([savePath 'Proj51_HPF.tif'],proj51_hpf);