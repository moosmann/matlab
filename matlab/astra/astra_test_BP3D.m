vol_geom = astra_create_vol_geom(128, 128, 128);

angles = linspace2(0, pi, 180);
proj_geom = astra_create_proj_geom('parallel3d', 1.0, 1.0, 128, 192, angles);

% Create a simple hollow cube phantom
cube = zeros(128,128,128);
cube(17:112,17:112,17:112) = 1;
cube(33:96,33:96,33:96) = 0;

% Create projection data from this
[~, proj_data] = astra_create_sino3d_cuda(cube, proj_geom, vol_geom);

% Display a single projection image
%figure, imshow(squeeze(proj_data(:,20,:))',[])

% Create a data object for the reconstruction
rec_id = astra_mex_data3d('create', '-vol', vol_geom);

% Set up the parameters for a reconstruction algorithm using the GPU
cfg = astra_struct('BP3D_CUDA');
cfg.ReconstructionDataId = rec_id;

%proj_geom = astra_create_proj_geom('parallel3d', 1.0, 1.0, 128, 192, pi/2);
% Free

% We generate the same geometry as the circular one above. 
angles = linspace2(0, pi, 1) + 0*pi/4;
disp( angles(1)*180/pi)
vectors = zeros(numel(angles), 12);
for i = 1:numel(angles)
  % ray direction
  vectors(i,1) = sin(angles(i));
  vectors(i,2) = -cos(angles(i));
  vectors(i,3) = 0;

  % center of detector
  %vectors(i,4:6) = 0;
  vectors(i,4) = cos( angles(i) );
  vectors(i,5) = sin( angles(i) );
  vectors(i,6) = 0;

  % vector from detector pixel (0,0) to (0,1)
  vectors(i,7) = cos(angles(i));
  vectors(i,8) = sin(angles(i));
  vectors(i,9) = 0;

  % vector from detector pixel (0,0) to (1,0)
  vectors(i,10) = 0;
  vectors(i,11) = 0;
  vectors(i,12) = 1;
end

% Parameters: #rows, #columns, vectors
proj_geom = astra_create_proj_geom('parallel3d_vec', 128, 192, vectors);

[proj_id, ~] = astra_create_sino3d_cuda(cube, proj_geom, vol_geom);
%cfg.ProjectionDataId = proj_id;
sino = proj_data(:,1,:);
%sino = permute( proj_data(:,1,:), [1 3 2]);
sino_id = astra_mex_data3d('create', '-proj3d', proj_geom, sino);
cfg.ProjectionDataId = sino_id;

% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run 150 iterations of the algorithm
% Note that this requires about 750MB of GPU memory, and has a runtime
% in the order of 10 seconds.
astra_mex_algorithm('iterate', alg_id, 150);

% Get the result
rec = astra_mex_data3d('get', rec_id);
f = figure('Name', sprintf( 'angles: %f', angles(1)*180/pi));
imagesc(squeeze(rec(:,:,65)));


% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data3d('delete', rec_id);
astra_mex_data3d('delete', proj_id);
