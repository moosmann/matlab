
%loadPath = '/home/astra05/katiN/Walnut3D/';
loadPath = '/home/jmoosmann/data/walnut/';
%loadName = 'Walnut_120_projections_384x328.mat'
loadName = 'Data82.mat';
load([loadPath loadName])


sinogram3D = reshape(sinogram3D,[384,328,120]);
sinogram3D = permute(sinogram3D,[2,3,1]);
[ndet_y,n_angles,ndet_x] = size(sinogram3D)

detSpace_x = 1.0;
detSpace_y = 1.0;
source_origin = 110; %mm 
origin_det = 300; %mm
angles = linspace(0, 2*pi, n_angles);
N = 128;
geomType = 'cone';
scale_y = ndet_y/114.8; % pix/mm
scale_x = ndet_x/115.2;


vol_geom = astra_create_vol_geom(N,N,N);
proj_geom = astra_create_proj_geom('cone', detSpace_x, detSpace_y, ndet_x, ndet_y, angles,source_origin*scale_x,origin_det*scale_x);


A = opTomo('cuda',proj_geom,vol_geom);
b = reshape(sinogram3D,[prod(size(sinogram3D)) 1]);

reco = lsqr(A,b,1e-4,1);
reco = reshape(reco,[vol_geom.GridRowCount vol_geom.GridColCount vol_geom.GridSliceCount]);
figure(1),clf,imagesc(reco(:,:,N/2)),colormap(gray)

B = opKron(opWavelet2(N,N,'Haar'),opWavelet(N,'Haar'));
reco_wave = lsqr(A*B',b,1e-4,2);
reco_wave = B'*reco_wave;
reco_wave = reshape(reco_wave,[vol_geom.GridRowCount vol_geom.GridColCount vol_geom.GridSliceCount]);

figure(2),clf,imagesc(reco_wave(:,:,N/2)),colormap(gray)
