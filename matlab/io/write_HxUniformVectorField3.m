function write_HxUniformVectorField3( filename, array, dtype, bounding_box, verbose )
% Write 3D volumetric vector field to Amira/Avizo compliant
% HxUniformVectorField3.
%
% ARGUMENTS
% filename : string. Output filename.
% array : 4D vector array, the first index defines the vector and the
%   latter three the shape of 3D volume in space.
% dtype : string, default: as input. data type. 'single' or 'double'
% bounding_box : 6 component vector, defines the size of the volume in
%   space. Default: [0, shape(1) - 1, 0, shape(2) - 1, 0, shape(3) - 1]
%
%
% Written by J. Moosmann, 2018-08-24. Last version: 2018-08-24

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    filename = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/HxUniformVectorField3/vec.am';
end
if nargin < 2
    array = ones( [3 100 100 100], 'single' );
end
if nargin < 3
    dtype = '';
end
if nargin < 4
    bounding_box = [];
end
if nargin < 5
    verbose = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = toc;

PrintVerbose( verbose, '\nWrite 4D vector array to Avizo/Amira-compliant HxUniformVectorField3' )

% Datatype
if isempty( dtype )
    dtype = class( array );
end

% Volumetric shape
shape = size( array );
shape(1) = [];

% Bounding box
if isempty( bounding_box )
    bounding_box = [0, shape(1) - 1, 0, shape(2) - 1, 0, shape(3) - 1];
end

%% Header
am_header = sprintf( [ ...
    sprintf('# Avizo BINARY-LITTLE-ENDIAN 2.1\n')...
    newline...
    newline...
    sprintf( 'define Lattice %u %u %u\n', shape)...
    newline...
    sprintf('Parameters {\n')...
    sprintf( '    Content "%ux%ux%u float[3], uniform coordinates",\n', shape)...
    sprintf( '    BoundingBox %u %u %u %u %u %u,\n', bounding_box)...
    sprintf( '    CoordType "uniform"\n' )...
    sprintf('}\n')...
    newline...
    sprintf('Lattice { float[3] Data } @1\n')...
    newline...
    sprintf('# Data section follows\n')...
    sprintf('@1\n')...
    ]);

%% Write Avizo/Amira format

% Open file stream
fid = fopen( filename, 'w' );

% Write header
fprintf( fid, '%s', am_header);

% Write array
fwrite( fid, array(:), dtype, 'l' );

% Close file
fclose( fid );

PrintVerbose( verbose, ', done in %.1f s = %.2f min', toc - t, (toc - t) / 60 )

% END
PrintVerbose( verbose, newline)

