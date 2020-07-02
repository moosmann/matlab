num = 1100;
bin = 4;

%% FULL
% intensity
att_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/processed/phase_1000/reco';
att = imread( sprintf( '%s/reco_%06u.tif', att_path, num ) );
attf = FilterHisto(Binning( att, bin ), 4);

% tie
tie_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/processed/phase_1000/reco_phase_tie_regPar2p50';
tie = imread( sprintf( '%s/reco_%06u.tif', tie_path, num ) );
tief = Binning( tie, bin);
tief = FilterHisto( tief, 4, 0.25);

% qp
qp_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/processed/phase_1000/reco_phase_qpcut_regPar2p50_binFilt0p200_cutoff2p00pi';
qp_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/scratch_cc/c20160803_001_pc_test/phase_1000/reco_phase/qpcut_regPar2p50_binFilt0p150_cutoff1p00pi';
qp = imread( sprintf( '%s/reco_%06u.tif', qp_path, num ) );
qpf = Binning( qp, bin);
qpf = FilterHisto( qpf, 4, 0.25);

%% ROI
x = IndexParameterToRange([0.4 0.6], size( att, 1));
y = IndexParameterToRange([0.4 0.6], size( att, 2));
nstd = 3;5;
roi = 0.25;
attr = FilterHisto( att(x,y), nstd, roi);
tier = FilterHisto( tie(x,y), nstd, roi);
qpr = FilterHisto( qp(x,y), nstd, roi);

mask = boolean(zeros( size( attf ) ));
xf = IndexParameterToRange([0.4 0.6], size( attf, 1));
yf = IndexParameterToRange([0.4 0.6], size( attf, 2));
mask( xf, [yf(1) - (1:4), yf(end) + (1:4)] ) = 1;
mask( [xf(1) - (1:4), xf(end) + (1:4)], yf ) = 1;

attf( mask ) = max( attf(:) );
tief( mask ) = max( tief(:) );
qpf( mask ) = max( qpf(:) );


%% Plot
n=2;
m=3;
h2 = figure(2);

subplot(n,m,1)
imsc( attf );
axis equal tight
title(sprintf('intensity'))
colorbar

subplot(n,m,2)
imsc( tief );
axis equal tight
title(sprintf('tie 2.5'))
colorbar

subplot(n,m,3)
imsc( qpf );
axis equal tight
title(sprintf('qp 2.5'))
colorbar

%% ROI
subplot(n,m,4)
imsc( attr );
axis equal tight
title(sprintf('intensity'))
colorbar

subplot(n,m,5)
imsc( tier );
axis equal tight
title(sprintf('tie 2.5'))
colorbar

subplot(n,m,6)
imsc( qpr );
axis equal tight
title(sprintf('qp 2.5'))
colorbar

%% save
%im_path = sprintf('/home/moosmanj/images/%s/%s/reco', beamtime_id, scan_name);
im_path = '/home/moosmanj/images/c20160803_001_pc_test/phase_1000/reco';

im = attf;
name = 'reco_int';
write32bitTIF(sprintf('%s/%s.tif', im_path, name), im);
imwritesc( im, sprintf('%s/%s.png', im_path, name))

im = tief;
name = 'reco_tie';
write32bitTIF(sprintf('%s/%s.tif', im_path, name), im);
imwritesc( im, sprintf('%s/%s.png', im_path, name))

im = qpf;
name = 'reco_qp';
write32bitTIF(sprintf('%s/%s.tif', im_path, name), im);
imwritesc( im, sprintf('%s/%s.png', im_path, name))

im = attr;
name = 'reco_int_roi';
write32bitTIF(sprintf('%s/%s.tif', im_path, name), im);
imwritesc( im, sprintf('%s/%s.png', im_path, name))

im = tier;
name = 'reco_tie_roi';
write32bitTIF(sprintf('%s/%s.tif', im_path, name), im);
imwritesc( im, sprintf('%s/%s.png', im_path, name))

im = qpr;
name = 'reco_qp_roi';
write32bitTIF(sprintf('%s/%s.tif', im_path, name), im);
imwritesc( im, sprintf('%s/%s.png', im_path, name))

