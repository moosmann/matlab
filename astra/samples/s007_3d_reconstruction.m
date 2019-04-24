% -----------------------------------------------------------------------
% This file is part of the ASTRA Toolbox
%
% Copyright: 2010-2018, imec Vision Lab, University of Antwerp
%            2014-2018, CWI, Amsterdam
% License: Open Source under GPLv3
% Contact: astra@astra-toolbox.com
% Website: http://www.astra-toolbox.com/
% -----------------------------------------------------------------------

clear all
aclear
% BP3D_CUDA
astra_mex('set_gpu_index', [0 1]);

vol_shape = [128 129 130]; sino_shape = [190 4097 200];
%vol_shape = [2047 2047 79]; sino_shape = [1280 6742 50];
sino_hor = sino_shape(1);
sino_proj = sino_shape(2);
sino_ver = sino_shape(3);

vol_geom = astra_create_vol_geom( vol_shape(2), vol_shape(1), vol_shape(3) );

angles = linspace2(0, pi, sino_proj);
proj_geom = astra_create_proj_geom('parallel3d', 1.0, 1.0, sino_ver, sino_hor, angles);

% Create a simple hollow cube phantom
cube = zeros( vol_shape(1), vol_shape(2), vol_shape(3) );
r1 = 0.15;
r2 = 1 -r1;
cube( ceil(r1*vol_shape(1)):floor(r2*vol_shape(1)),ceil(r1*vol_shape(2)):floor(r2*vol_shape(2)), ceil(r1*vol_shape(3)):floor(r2*vol_shape(3))) = 1;

r1 = 0.25;
r2 = 1 -r1;
cube( ceil(r1*vol_shape(1)):floor(r2*vol_shape(1)),ceil(r1*vol_shape(2)):floor(r2*vol_shape(2)), ceil(r1*vol_shape(3)):floor(r2*vol_shape(3))) = 0;

% Create projection data from this
[proj_id, proj_data] = astra_create_sino3d_cuda(cube, proj_geom, vol_geom);

% Display a single projection image
figure(1), imshow(squeeze(proj_data(:,20,:))',[])

% Create a data object for the reconstruction
rec_id = astra_mex_data3d('create', '-vol', vol_geom);

%cfg = astra_struct('SIRT3D_CUDA');

%% BP3D
cfg = astra_struct('BP3D_CUDA');
fprintf( '\n%s', cfg.type )
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id, 10);

%% FP3D
if 0
    cfg = astra_struct('FP3D_CUDA');
    fprintf( '\n%s', cfg.type )
    cfg.VolumeDataId = rec_id;
    cfg.ProjectionDataId = proj_id;
    alg_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('iterate', alg_id, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the result
rec = astra_mex_data3d('get', rec_id);
figure(2), imshow(squeeze(rec(:,:,65)),[]);

ainfo

% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data3d('delete', rec_id);
astra_mex_data3d('delete', proj_id);
