% 
sino=Readstack('/export/scratch1/moosmann/art/Batenburg/MediumQualityDataSet/sino',100);
% Parameters from PyHST par file needed for tomographic reconstruction
%NUM_FIRST_IMAGE = 1 # No. of first projection file
%NUM_LAST_IMAGE  = 832 # No. of last projection file
%NUM_IMAGE_1 = 1488 # Number of pixels horizontally
%NUM_IMAGE_2 = 1300 # Number of pixels vertically
%ANGLE_BETWEEN_PROJECTIONS = 0.21635 # Increment angle in degrees
%ROTATION_AXIS_POSITION = 725.500000 # Position in pixels
rotpos = 725.5;
%END_VOXEL_1 = 1488 # X-end of reconstruction volume
%END_VOXEL_2 = 1488 # Y-end of reconstruction volume
%END_VOXEL_3 = 1300 # Z-end of reconstruction volume
sino = sino(:,1:round(2*rotpos));
sino = 100*(sino-mean(sino(:)));
%sino = (sino(1:12:end,round(size(sino,2)/2)+(-256:255)));
%sino = (sino(1:1:end,round(size(sino,2)/2)+(-512:511)));


vol_geom = astra_create_vol_geom([size(sino,2),size(sino,2)]);
proj_geom = astra_create_proj_geom('parallel',1.0,size(sino,2),linspace2(0,pi,size(sino,1)));

tic
[v_fbp_id, v_fbp]=astra_create_reconstruction_cuda('SIRT_CUDA',proj_geom,vol_geom,sino,1,0,0,0,0,0,0);
proj_id = astra_create_projector('line', proj_geom, vol_geom);
%[v_fbp_id, v_fpb]=astra_create_reconstruction('FBP',proj_id,sino,1,0,0,0,0,0,0);
%[v_fbp_id, v_fbp] = astra_create_fbp_reconstruction(sino, proj_id);
toc
dimx = size(v_fbp,1);
x = round(dimx/2)+(-ceil(dimx/2/sqrt(2)):floor(dimx/2/sqrt(2)));
itool(v_fbp(x,x))