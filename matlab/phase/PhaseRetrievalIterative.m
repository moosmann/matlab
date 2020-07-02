proj_path = '/asap3/petra3/gpfs/p07/2019/data/11006991/processed/syn008_Ti_12w_47R_z0030_t135ms/flat_corrected/rawBin4';

tic
fprintf( '\n Create datastore ' )
ds = datastore( proj_path, ...
    'IncludeSubfolders', true, ...
    'FileExtensions', '.tif', ...
    'Type', 'image');
fprintf( ' in %g s', toc )
s = dir( [proj_path '/*tif' ] );
toc
