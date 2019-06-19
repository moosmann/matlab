
%% Reference images DCM CCD
% folder = '/asap3/petra3/gpfs/p05/2018/data/11004263/raw/syn004_96R_Mg5Gd_8w_a';
% pattern = '*ref*';
% step_size = 1;
% rect2 = [705 14 1247 3036];
% vol_transform = 'rot90( v )';
% outpath = '/asap3/petra3/gpfs/common/p05/jm/pictures';
% outname = 'proj_dcm_ccd';
% res_xy = 1020;
% res_z = 64;

% %% Projections
% folder = '/asap3/petra3/gpfs/p05/2018/data/11004263/raw/syn004_96R_Mg5Gd_8w_a';
% pattern = '*img*';
% step_size = 1:100;
% rect2 = [705 14 1247 3036];
% vol_transform = 'rot90( v )';
% outpath = '/asap3/petra3/gpfs/common/p05/jm/pictures';
% outname = 'proj_dcm_ccd';
% res_xy = 1020;
% res_z = 64;

%% Reference images
folder = '/asap3/petra3/gpfs/p05/2018/data/11004936/raw/syn005_55R_Mg5Gd_12w_load_000';
pattern = '*ref*';
step_size = 1;
rect2 = [];
vol_transform = 'rot90( v )';
outpath = '/asap3/petra3/gpfs/common/p05/jm/pictures';
outname = 'proj_dmm_ccd';
res_xy = 720;
res_z = 64;

% %% Phase retrieval rot axis
% outpath = '/asap3/petra3/gpfs/common/p05/jm/pictures';
% outname = 'wood_find_rot_axis';
% vol_transform = 'rot90( v )';
% rect2 = [1457         857         360        1383];
% res_xy = 1020;
% res_z = 64;
% 
% %% Phase retrieval reg par
% outpath = '/asap3/petra3/gpfs/common/p05/jm/pictures';
% outname = 'wood_find_reg_par';
% vol_transform = '';
% rect2 = [];
% res_xy = 1020;
% res_z = 64;

tic
%% Read
if ~exist( 'vol', 'var' )
    vol = read_images_to_stack( folder, step_size, pattern );
end

%% Crop
if isempty( rect2 )
    im = double( sum( vol, 3 ) );
    im = FilterOutlier( im, 0.02 );
    im = normat( im );
    figure( 'WindowState', 'maximized' )
    [imc, rect2] = imcrop( im );
    rect2 = round( rect2 );
    fprintf( '\n Cropping rectangle: \n' )
    fprintf( ' rect2 = [ %u %u %u %u ];\n', rect2 )
end
xx = rect2(1) + (0:rect2(3) -1);
yy = rect2(2) + (0:rect2(4) -1);

volc = vol( yy, xx, : );

%% Transform
if ~isempty( vol_transform )
    fprintf( '\nTransforming' )
    v = volc;
    volc = eval( vol_transform );
end

%% Rescale
fprintf( '\nRescaling/Binning: ' )
t = toc;
[x, y, z] = size( volc );
dx = max( round( max( x, y) / res_xy ), 1 );
if length( res_z ) == 1
    dz = max( floor( z / res_z ), 1 );
    zz = 1:dz:z;
else
    zz = res_z;
end
volr = Binning( double( volc(:,:,zz) ), [dx dx 1] );
% dy = dx;
% dz = 1;
% Xq = 1:dx:x;
% Yq = 1:dy:y;
% Zq = 1:dz:z;
% volr = interp3( vol, Xq, Yq, Zq, 'nearest' );
fprintf( '%.0f s', toc -t )
fprintf( '\nVolume size: %u %u %u', size( volr ) )

%% Adjust dynamic
%vola = FilterOutlier( volr, 0.02 );


%% Write
fprintf( '\nWrite GIF' )
if isempty( outpath )
    outpath = regexprep( folder, 'raw', 'scratch_cc' );
end
CheckAndMakePath( outpath )
if isempty( outname )
    outname = regexprep( pattern, '*', '' );
end
filename = sprintf( '%s/%s.gif', outpath, outname );
fprintf( '\noutput file path:\n %s', filename )

write_gif( volr, filename )
% End
fprintf( '\nFINISHED in %.0f s \n', toc )
