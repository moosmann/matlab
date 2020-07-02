function write_divergence_vector_tiffs_to_HxUniformScalarField( input_folder, steps_pattern, vector_pattern, out_path, out_name, steps_to_process )
% Read in vector data stored as single component 3D volumes, compute
% divergence and save as Amira/Avizo HxUniformScalarField
%
% Written by J. Moosmann, 2018-08-23. Last version: 2018-08-30

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/3d_tiff';
end
if nargin < 2
    steps_pattern = 'stk_*_x.tif';
end
if nargin < 3
    vector_pattern = {'vec_x', 'vec_y', 'vec_z'};
end
if nargin < 4
    out_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/HxUniformScalarField';
end
if nargin < 5
    out_name = 'div';
end
if nargin < 6
    steps_to_process = [];
end
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( '\nRead single component vector tiffs, compute divergence and save as Avizo/Amira-compliant HxUniformVectorField3' )
fprintf( '\n output folder: %s', out_path)
CheckAndMakePath( out_path );
CheckTrailingSlash( out_path );

CheckTrailingSlash( input_folder );

fs = dir( [input_folder steps_pattern] );
num_data_sets = numel( fs );
fprintf( '\n number of data sets found: %u', num_data_sets )

if isempty( steps_to_process )
    steps_to_process = 1:num_data_sets;
end

filename = [fs(1).folder filesep fs(1).name];
tiff_info = imfinfo( filename );

shape = [ tiff_info(1).Height tiff_info(1).Width numel( tiff_info ) ];

vec = zeros( [ 3 shape], 'single' );
for nn = 1:numel( steps_to_process )
    
    % Read data
    
    t = toc;
    num = steps_to_process(nn);
    fprintf( '\nStep: %u. Reading data', num )
    
    filename = [input_folder fs(nn).name ];
    vec(1,:,:,:) = read_multitif( filename, [], [], tiff_info, 0 );
    
    filename = regexprep( filename, vector_pattern{1}, vector_pattern{2} );
    vec(2,:,:,:) = read_multitif( filename, [], [], tiff_info, 0 );
    
    filename = regexprep( filename, vector_pattern{2}, vector_pattern{3} );    
    vec(3,:,:,:) = read_multitif( filename, [], [], tiff_info, 0 );
    
    fprintf( ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )
    
    % Divergence
    fprintf( '\n         Compute divergence, step: %u', num )
    div = divergence( vec );
    fprintf( ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )

    % Write    
    t = toc;
    fprintf( '\n         Writing data, step: %u', num )
    
    % Write Amira/Avizo binary
    name = sprintf( '%s_%02u.am', out_name, num );
    filename = sprintf( '%s%s', out_path, name );
    fprintf( ', output filename: %s', name)
    dtype = class( vec );
    write_HxUniformScalarField( filename, div, dtype, [], 0 )
            
    fprintf( ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )
    
end
fprintf( '\nTotal processing time %.1f s = %.2f min', toc, toc  / 60 )


% END
fprintf( '\n' )
