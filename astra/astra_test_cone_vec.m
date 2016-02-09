% Create a 3D volume geometry.
% Parameter order: rows, colums, slices (y, x, z)
vol_geom = astra_create_vol_geom(50, 50 , 50);

% initialized to 3.0
data = astra_mex_data3d('create', '-vol', vol_geom, 3.0);


% Projection geometry
det_row_count = 100;
det_col_count = 100;

x = 0.70710678;
vectors = [ 1 0 0  0 0 0   0 1 0  0 0 1 ; 
            x x 0  0 0 0  -x x 0  0 0 1;
            0 1 0  0 0 0  -1 0 0  0 0 1 ];
new_ind = [2, 1, 0, 5, 4, 3, 8, 7, 6, 11, 10, 9] + 1;
v = zeros(size(vectors));
for nn = 1:size(vectors,1)
    v(nn, :) = vectors(nn, new_ind);
end

proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, vectors);
proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, v);


[proj_id, proj_data] = astra_create_sino3d_cuda(data, proj_geom, vol_geom);


% Retrieve data:
%R = astra_mex_data3d('get', v1);

% Retrieve data as a single array. Since astra internally stores
% data as single precision, this is more efficient:
%R = astra_mex_data3d('get_single', v1);



% Delete all created data objects
%astra_clear

ishow(squeeze(proj_data(:,1,:)))
ishow(squeeze(proj_data(:,2,:)))
ishow(squeeze(proj_data(:,3,:)))