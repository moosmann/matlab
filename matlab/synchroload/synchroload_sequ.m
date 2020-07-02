clear all
tic
t = toc;
p = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load_';
% t0: 722
% t late: 653
outpath = '/home/moosmanj/images/load_sequ';

sln = 653;

ind = [0:2:18];


for nn = length( ind ):-1:1
    num = ind(nn);
    filename = sprintf( '%s%02u/flat_corrected/proj_%06u.tif', p, num, sln);    
    %fprintf( '\n%u : %s', nn, filename )
    proj(:,:,nn) = imread( filename, 'tif' );
    %read_tif( filename );
end

for nn = length( ind ):-1:2
    out = ImageCorrelation( proj(:,:,nn-1),proj(:,:,nn), 0 );
    x(nn) = out.shift1;
    y(nn) = out.shift2;    
end

xx = cumsum( x );
yy = cumsum( y );


for nn = length( ind ):-1:1
    num = ind(nn);
    filename = sprintf( '%s%02u/reco/float_binned/reco_%06u.tif', p, num, sln + 0*round(xx(nn)));    
    fprintf( '\n%u : %s', nn, filename )
    v(:,:,nn) = imread( filename, 'tif' );
    %read_tif( filename );
end

%% Write
filename = sprintf( '%s/load_sequ_z%04u.tif', outpath, sln );
imwrite( v, filename, 'tif' )

fprintf( '\n Elapsed time : %.0f s (%.1f min)', toc - t, (toc - t)/60 )

%% Show
volshow( v, 3, 1 )