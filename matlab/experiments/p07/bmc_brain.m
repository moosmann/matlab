%% phase retrieval 2D
inpath = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/trans02_180';
sino = read_images_to_stack(inpath,1,'*.tif',[],1,1);

 EnergyDistancePixelsize = [67e3 0.8 4*1.27e-6];



%% phase retrieval 3D
%% Read data
fprintf('\nReading volume')
%fn = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco/noNegLog_ringFiltJM/float_rawBin2/bmc006_B2_8rings_h2.h5';
%fn = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco/noNegLog/float_rawBin2/bmc006_B2_8rings_h2_dering_rw040.h5';%
fn = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco/noNegLog/float_rawBin2/bmc006_B2_8rings_h2.h5';
fn = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco/noNegLog_ringFiltJM/float_rawBin2/bmc006_B2_8rings_h2_dering_rw100.h5';
fn = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco/noNegLog/float_rawBin2/bmc006_B2_8rings_h2.h5';
[~,scan_name] = fileparts(fn);
fprintf('\n input volume filename:\n %s',fn)
tic;
vol = h5read(fn,'/volume');
t = toc;
fprintf(' done in %.0f s',t)

%% Phase Retrieval
fprintf('\nPhase retrieval')
pha = -PhaseRetrieval3D(vol);

%% Write data
fprintf('\nWriting volume')
tic;
%outpath = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco_rf_phase3D/float_rawBin2/';
outpath = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco_phase3D/float_rawBin2/';
outpath = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco_rfjm_rr_phase3D/float_rawBin2/';
outpath = '/asap3/petra3/gpfs/p07/2023/data/11016192/processed/bmc006_B2_8rings_h2/reco_phase3D/float_rawBin2/';
fprintf('\n output folder: /n%s', outpath)
CheckAndMakePath(outpath);
fn = sprintf('%s%s.h5',outpath,scan_name);
fprintf('\n output file name:\n %s',fn)
h5create(fn,'/volume',size(pha),'Datatype','single')
h5write(fn,'/volume',pha)
t = toc;
fprintf('\n done in %.0f s',t)
fprintf('\n')

