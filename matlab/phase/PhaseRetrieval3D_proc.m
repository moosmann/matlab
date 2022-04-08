reco_path = pwd;

fprintf( '\nReading volume' )
if ~exist('vol','var')
    vol = read_images_to_stack(reco_path,1,'*tif',[],1,1);
end

fprintf( '\nBinning' )
vol2 = Binning(vol,2);

method = 'tie';
EnergyDistancePixelsize = [57e3 1.4 2*5*0.9e-6];
reg_par = 2;
bin_filt = 0.1;
padding = 0;

fprintf( '\nPhase retrieval 3D' )
pha1 = PhaseRetrieval3D(vol, method, EnergyDistancePixelsize, reg_par,bin_filt,padding);

pha2 = PhaseRetrieval3D(vol2, method, EnergyDistancePixelsize, reg_par,bin_filt,padding);

outpath = '/asap3/petra3/gpfs/p07/2022/data/11014410/processed/hereon04_mhh_ts02/phase3d/tie_regPar2p0/float_rawBin5_noPadding';
CheckAndMakePath(outpath)
WriteVol(pha1,outpath,'hereon04_mhh_ts02_',3);