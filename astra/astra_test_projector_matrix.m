clear all
vol_geom = astra_create_vol_geom(256, 256);
proj_geom = astra_create_proj_geom('parallel', 1.0, 20, linspace(0,pi,18));

proj_id = astra_create_projector('line', proj_geom, vol_geom);
proj_id_cu = astra_create_projector('cuda', proj_geom, vol_geom);

mat_id = astra_mex_projector('matrix', proj_id);
mat_id_cu = astra_mex_projector('matrix', proj_id_cu);

fprintf('\nproj_id is cuda: %u', astra_mex_projector('is_cuda', proj_id))
fprintf('\nproj_id is cuda: %u', astra_mex_projector('is_cuda', proj_id_cu))

smat = astra_mex_matrix('get', mat_id);
smat_cu = astra_mex_matrix('get', mat_id_cu);

mat = full(smat);
mat_cu = full(smat_cu);

P = phantom(256)';
s = smat * P(:);
s_cu = smat_cu * P(:);

s = reshape(s, [proj_geom.DetectorCount size(proj_geom.ProjectionAngles, 2)])';
s_cu = reshape(s_cu, [proj_geom.DetectorCount size(proj_geom.ProjectionAngles, 2)])';

ishow(s)
ishow(s_cu)

% ainfo
astra_mex_matrix('delete', mat_id)
astra_mex_matrix('delete', mat_id_cu)


fprintf('\n')
