

InPath = '/mnt/LSDF/users/moosmann/CWI_DATA/ESRF_MI1079_ID19_July2011_inlineTomo/vol/Xenopus_4cell_20keV/sino_qpRP25BF01Proj1__1600x1962x0100';
OutPath = '/mnt/LSDF/users/moosmann/phd/figures/TomoPadding/';
OutPathTif = '/mnt/LSDF/users/moosmann/phd/figures/TomoPadding/tif/';
CheckAndMakePath(OutPathTif)

%% Read sino
sino0 = double( imread( [InPath '/slice_0050.tif'] ) );
sino  = SubtractMean(sino0);

%% Pad sino
sinopz = padarray(SubtractMean(sino),[0 981],0,'both');
sinopr = padarray(SubtractMean(sino),[0 981],'replicate','both');
sinops = padarray(SubtractMean(sino),[0 981],'symmetric','both');

%% Reco
vol0   = astra_make_fbp(sino0,2*pi*(0:1599)/1600,'FBP_CUDA');
vol   = astra_make_fbp(SubtractMean(sino),2*pi*(0:1599)/1600,'FBP_CUDA');
volpz = astra_make_fbp(SubtractMean(sinopz),2*pi*(0:1599)/1600,'FBP_CUDA');
volpr = astra_make_fbp(SubtractMean(sinopr),2*pi*(0:1599)/1600,'FBP_CUDA');
volps = astra_make_fbp(SubtractMean(sinops),2*pi*(0:1599)/1600,'FBP_CUDA');

%% Crop
x = 981+(1:1962);
volpz = volpz(x,x);
volpr = volpr(x,x);
volps = volps(x,x);

%% Save Tiff
write32bitTIF( [ OutPathTif 'tomo_qpRP25BF01_sinoPadNo_fbp.tif']    , vol );
write32bitTIF( [ OutPathTif 'tomo_qpRP25BF01_sinoPadZero_fbp.tif' ] , volpz );
write32bitTIF( [ OutPathTif  'tomo_qpRP25BF01_sinoPadRep_fbp.tif' ] , volpr );
write32bitTIF( [ OutPathTif 'tomo_qpRP25BF01_sinoPadSym_fbp.tif' ]  , volps );

%% Save PNG
WriteImage( [ OutPath 'tomo_qpRP25BF01_sinoPadNo_fbp' ]   , Binning(vol)  , 'png' );
WriteImage( [ OutPath 'tomo_qpRP25BF01_sinoPadZero_fbp' ] , Binning(volpz) , 'png' );
WriteImage( [ OutPath  'tomo_qpRP25BF01_sinoPadRep_fbp']  , Binning(volpr) , 'png' );
WriteImage( [ OutPath 'tomo_qpRP25BF01_sinoPadSym_fbp' ]  , Binning(volps) , 'png' );

%% Show
itool(vol0),itool(vol),itool(volpz),itool(volpr),itool(volps)
cd(OutPath)