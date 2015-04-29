%sino=Readstack('/export/scratch1/moosmann/art/Batenburg/LowQualityDataSet/sino_HigherCounts',100);
%rotpos = 626;%725.5;
filename = '/export/scratch1/moosmann/art/Batenburg/MediumQualityDataSet/sino_filtSinoTriangle/sino_0610.tif';
%filename = '/export/scratch1/moosmann/art/Batenburg/MediumQualityDataSet/sino/sino_0610.tif';
sino = imread(filename);
rotpos = 725.5000;

sino = sino(:,1:round(2*rotpos));
%sino = (sino(1:12:end,round(size(sino,2)/2)+(-256:255)));
%sino = (sino(1:1:end,round(size(sino,2)/2)+(-512:511)));
sino = sino(1:4:end,1:2:end);

vol_geom = astra_create_vol_geom([size(sino,2),size(sino,2)]);
proj_geom = astra_create_proj_geom('parallel',1.0,size(sino,2),linspace2(0,pi,size(sino,1)));


tic
[rec_id, rec]=astra_create_reconstruction_cuda('SIRT_CUDA',proj_geom,vol_geom,sino,200,0,0,0,0,0,0);

%proj_id = astra_create_projector('line', proj_geom, vol_geom);
%[v_fbp_id, v_fpb]=astra_create_reconstruction('FBP',proj_id,sino,1,0,0,0,0,0,0);
%[rec_id, rec] = astra_create_fbp_reconstruction(sino, proj_id);
toc

dimx = size(rec,1);
x = round(dimx/2)+(-ceil(dimx/2/sqrt(2)):floor(dimx/2/sqrt(2)));
itool(rec(x,x))