clear
aclear

% Volume shape
vol_shape = [10, 20, 30];
vol_size = [-1, 1, -2, 2, -3, 3];


% Volume shape and size
row_count = vol_shape(2);
col_count = vol_shape(1);
slice_count = vol_shape(3);
min_x = vol_size(3);
max_x = vol_size(4);
min_y = vol_size(1);
max_y = vol_size(2);
min_z = vol_size(5);
max_z = vol_size(6);

% Matlab array
vol = ones(vol_shape,'single');

% ASTRA array
%vol_shape_astra = [vol_shape(2), vol_shape(1), vol_shape(3)];
%vol_geom = astra_create_vol_geom(vol_shape_astra);
vol_geom = astra_create_vol_geom(row_count, col_count, slice_count, min_x, max_x, min_y, max_y, min_z, max_z);
vol_id_link = astra_mex_data3d_c('link', '-vol', vol_geom, vol, 1, 0);

% proj
detector_spacing_x = 0.9;
detector_spacing_y = 1.1;
det_row_count = 20;
det_col_count = 50;
angles = 60;
proj_geom = astra_create_proj_geom('parallel3d', detector_spacing_x, detector_spacing_y, det_row_count, det_col_count, angles);
sino_id = astra_mex_data3d('create', '-proj3d', proj_geom);

ainfo