%% 2020 Sept P05



p.regdir = 'x';
p.steps = [1 2 5 6 7 3 4];
p.register = 0;
p.outlier_thresh = 0.0001;
p.proc_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/processed/';
p.voxel_size = 5.2e-6;
p.barcol = 'white';

p.scan_name = 'phase';

p.reco_sub = 'reco/float_rawBin4';
p.reco_sub = 'reco_phase/tie_regPar2p00/float_rawBin4';
vol = p05_load_sequ( p );
