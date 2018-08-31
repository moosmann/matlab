function write_divergence_vector_binaries_to_HxUniformScalarField( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, out_path, out_name, steps_to_process )
% Read in vector data stored as single component 3D volumes, compute
% divergence and save as Amira/Avizo HxUniformScalarField
%
% Written by J. Moosmann, 2018-08-30. Last version: 2018-08-30

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/raw';
end
if nargin < 2
    folder_pattern = 'out*';
end
if nargin < 3
    filename_pattern = 'stk*.raw';
end
if nargin < 4
    shape = [700 600 500];
end
if nargin < 5
    dtype = 'single';
end
if nargin < 6
    machinefmt = 'b';
end
if nargin < 7
    out_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/AvizoTest/HxUniformScalarField';
end
if nargin < 8
    out_name = 'div';
end
if nargin < 9
    steps_to_process = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

fprintf( '\nRead single component binary vector data, compute divergence and save as Avizo/Amira-compliant HxUniformScalarField' )
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

for nn = 1:numel( steps_to_process )
    
    %% Read vector field
    t = toc;
    num = steps_to_process(nn);
    fprintf( '\nStep %u: Reading data', num )
    
    folder = [fs(num).folder filesep fs(num).name];
    fs2 = dir( [folder filesep filename_pattern] );
    
    % Read vector field component 1
    mm = 1;
    filename = [fs2(mm).folder filesep fs2(mm).name ];
    vx = read_binary_data( filename, shape, dtype, machinefmt );
    
    % Read vector field component 2
    mm = 2;
    filename = [fs2(mm).folder filesep fs2(mm).name ];
    vy = read_binary_data( filename, shape, dtype, machinefmt );
    
    % Read vector field component 3
    mm = 3;
    filename = [fs2(mm).folder filesep fs2(mm).name ];
    vz = read_binary_data( filename, shape, dtype, machinefmt );
    
    fprintf( ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )
    
    %% Divergence
    fprintf( '\n         Compute divergence' )
    div = divergence( vx, vy, vz );
    fprintf( ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )
    
    %% Write Avizo/Amira format        
    t = toc;
    fprintf( '\n         Writing data' )
    
    name = sprintf( '%s_%02u.am', out_name, num );
    filename = sprintf( '%s%s', out_path, name );
    fprintf( ', output filename: %s', name)
    write_HxUniformScalarField( filename, div, dtype, [], 0 )
    
    fprintf( ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )
    
end
fprintf( '\nTotal processing time %.1f s = %.2f min', toc, toc  / 60 )

% END
fprintf( '\n' )
