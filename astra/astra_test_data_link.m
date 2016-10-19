clear
aclear

vol_shape = [3,4,2];
vol = ones(vol_shape,'single');

% vol
vol_shape_astra = [vol_shape(2), vol_shape(1), vol_shape(3)];
vol_geom = astra_create_vol_geom(vol_shape_astra);
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