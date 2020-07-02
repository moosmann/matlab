function write_single_component_vector_binaries_to_HxUniformVectorField3( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, out_path, out_name, steps_to_process )
% Read in vector data stored as single component 3D volumes and save as
% Amira/Avizo HxUniformVectorField3
%
%
% Written by J. Moosmann, 2018-08-23. Last version: 2018-08-24

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc';
end
if nargin < 2
    folder_pattern = 'out*';
end
if nargin < 3
    filename_pattern = 'vec*.vol';
end
if nargin < 4
    shape = [1240 1120 20];
end
if nargin < 5
    dtype = 'single';
end
if nargin < 6
    machinefmt = 'l';
end
if nargin < 7
    out_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc/HxUniformVectorField3';
end
if nargin < 8
    out_name = 'vec';
end
if nargin < 9
    steps_to_process = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

fprintf( '\nRead single component binary vector data and save as Avizo/Amira-compliant HxUniformVectorField3' )
fprintf( '\n output folder: %s', out_path)
CheckAndMakePath( out_path );
CheckTrailingSlash( out_path );

%% Read data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CheckTrailingSlash( input_folder );

fs = dir( [input_folder folder_pattern] );
num_data_sets = numel( fs );
fprintf( '\n number of data sets found: %u', num_data_sets )

if isempty( steps_to_process )
    steps_to_process = 1:num_data_sets;
end

v = zeros( [ 3 shape], 'single' );
for nn = 1:numel( steps_to_process )
    
    
    t = toc;
    num = steps_to_process(nn);
    fprintf( '\nStep: %u. Reading data', num )
    
    folder = [fs(num).folder filesep fs(num).name];
    fs2 = dir( [folder filesep filename_pattern] );
    
    % Read vector field component 1
    mm = 1;
    filename = [fs2(mm).folder filesep fs2(mm).name ];
    v(1,:,:,:) = read_binary_data( filename, shape, dtype, machinefmt );
    
    % Read vector field component 2
    mm = 2;
    filename = [fs2(mm).folder filesep fs2(mm).name ];
    v(2,:,:,:) = read_binary_data( filename, shape, dtype, machinefmt );
    
    % Read vector field component 3
    mm = 3;
    filename = [fs2(mm).folder filesep fs2(mm).name ];
    v(3,:,:,:) = read_binary_data( filename, shape, dtype, machinefmt );
    
    fprintf( ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )

    % Filter NAN
%     m = isnan(v);
%     num_nan = sum(m(:));
%     fprintf( ', NANs: %u (%f %%)', num_nan, num_nan / numel(m))
%     v(m) = 0;
%     m = isinf(v);
%     num_inf = sum(m(:));
%     fprintf( ', INFs: %u (%f %%)', num_inf, num_inf / numel(m))
%     v(m) = 0;
%     
    %% Write Avizo/Amira format    
    t = toc;
    fprintf( '\n         Writing data, step: %u', num )
    name = sprintf( '%s_%02u.am', out_name, num );
    filename = sprintf( '%s%s', out_path, name );
    fprintf( ', output filename: %s', name)
    write_HxUniformVectorField3( filename, v, dtype, [], 0 )
            
    fprintf( ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )
end
fprintf( '\nTotal processing time %.1f s = %.2f min', toc, toc  / 60 )


% END
fprintf( '\n' )
