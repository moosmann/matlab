clear all
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag' );

raw_path = '/asap3/petra3/gpfs/p05/2017/data/11002839/raw/ehh_2017_019_f/';
%raw_path = '/home/moosmanj/testTiff/';
filename = [raw_path 'proj_0000.tif'];

% Tiff info struct
tinfo = imfinfo( filename );
width = tinfo.Width;
height = tinfo.Height;

fprintf( '\n Width : %u', width )
fprintf( '\n Height : %u', height )
fprintf( '\n RowsPerStrip : %u', tinfo.RowsPerStrip )
fprintf( '\n PlanarConfiguration : %s', tinfo.PlanarConfiguration )

fprintf( '\n' )
num_proj = 10;
proj = zeros( height, width, num_proj, 'uint16');

meth = 'jm';
%meth = 'ri';
%meth = 'ml';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
for nn = 1:num_proj
    filename = sprintf( '%sproj_%04u.tif', raw_path, nn -1 );
    switch meth
        case 'jm'
            TifLink = Tiff( filename, 'r');
            %TifLink.setTag('Photometric', 1);
            proj(:,:,nn) = TifLink.read();
            TifLink.close();            
        case 'ml'
            proj(:,:,nn) = imread( filename, 'tif', 'Info', tinfo );
        case 'ri'
            proj(:,:,nn) = read_image( filename, 'tif' );
    end
end
t = toc;
fprintf( 'Elapsed time: %.1f s. Per image: %.3f\n', t, t / num_proj )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TifLink = Tiff( filename, 'r');
% %TifLink.Photometric = 1;
% strip_number = computeStrip( TifLink,  floor(height / 2 ));
% 
% fprintf( '\n ImageLength : %u', TifLink.getTag('ImageLength') )
% fprintf( '\n strip number : %u', strip_number )
% 
% s = TifLink.readEncodedStrip( strip_number );
% fprintf( '\n strip size : %u %u ', size( s ) )
% 
% TifLink.close();
