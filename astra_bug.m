function astra_bug(dimRange)

if nargin < 1
    dimRange = 1:50;
end

for dim = dimRange
    % Create volume geometry
    vol_geom = astra_create_vol_geom( [dim, dim, dim]);
    
    % Create volume
    rec_id = astra_mex_data3d('create', '-vol', vol_geom);
    
    % Intialize
    astra_mex_data3d('set', rec_id, 1 );
    
    % Fetch data object as 3D matrix
    % Results differ if data is fetched or not
    rec = astra_mex_data3d('get_single', rec_id);
    
    % clear astra memory
    astra_mex_data3d('delete', rec_id);
    astra_mex_data3d('clear');
    
    % Create volume geometry
    vol_geom = astra_create_vol_geom( [dim, dim, dim]);
    
    % Create volume
    rec_id = astra_mex_data3d('create', '-vol', vol_geom);
    
    % Fetch data object as 3D matrix
    rec = astra_mex_data3d('get_single', rec_id);
    
    % Check if volume is intialized to zero
    fprintf( '%g ', sum( rec(:) ) == 0 )
    
end
fprintf('\n')