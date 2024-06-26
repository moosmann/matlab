

if ~exist('vol','var')
    fprintf( '\nReading volume' )
    reco_path = pwd;
    vol = read_images_to_stack(reco_path,1,'*tif',[],1,1);
end

method = 'tie';
EnergyDistancePixelsize = [67e3 0.8 2*1.27e-6];%[57e3 1.4 2*5*0.9e-6];
reg_par = 1;
bin_filt = 0.1;
padding = 0;

fprintf( '\nPhase retrieval 3D' )
pha1 = PhaseRetrieval3D(vol, method, EnergyDistancePixelsize, reg_par,bin_filt,padding);

if 0
    fprintf( '\nBinning' )
    vol2 = Binning(vol,2);
    pha2 = PhaseRetrieval3D(vol2, method, EnergyDistancePixelsize, reg_par,bin_filt,padding);
end

fprintf('\nSaving tifs')
outpath = '/asap3/petra3/gpfs/p07/2023/data/11019133/scratch_cc/kulvait/wd_mouse_008/bmc008_85679_mousebrain_202311_DESY_height_g/recb2_ext_wlssd_tie1p0';
CheckAndMakePath(outpath)
WriteVol(pha1,outpath,'algvk_ext_avg_phase_',3); 