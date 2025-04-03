%clear proj

pp = gcp( 'nocreate' );

if isempty( pp )
    pp = parpool( 'Processes' );
end
num_workers = pp.NumWorkers;
fprintf( '\n\nWRITE/READ TEST' )
if ~exist( 'nn', 'var' )
    nn = 0;
end

fprintf( '\n parpool workers %u', num_workers )

out_path = '/data/hereon/wp/user/moosmanj/test_io';
CheckAndMakePath( out_path );
CheckAndMakePath( [ out_path '/hdf'] )
CheckAndMakePath( [ out_path '/tif_multi'] )
CheckAndMakePath( [ out_path '/tif_single'] )
CheckAndMakePath( [ out_path '/binary'] )

%raw_path = '/asap3/petra3/gpfs/p07/2020/data/11010172/raw/swerim_21_12_oh_a/tiff0000/';
raw_path = '/asap3/petra3/gpfs/p07/2022/data/11012618/processed/etna_068_cse100222g_f/reco/float_rawBin2/';
fs = dir( [raw_path '*tif'] );
fprintf( '\n %u tiff files found in\n  %s\n', numel( fs ), raw_path )

fn = [fs(1).folder filesep fs(1).name];
imf = imfinfo( fn, 'tif' );
imsize = [imf.Height imf.Width];

num_proj = numel( fs );
if ~exist( 'proj', 'var' )
    fprintf( '\nReading images' )
    tic
    %proj = zeros( [imsize num_proj], 'uint16' );
    proj = zeros( [imsize num_proj], 'single' );
    parfor n = 1:num_proj
        fn = [fs(n).folder filesep fs(n).name];
        proj(:,:,n) = imread( fn, 'tif', 'Info', imf )
    end
    t = toc;
end
ishow( proj(:,:,round(size(proj,3)/2)))
ishow( proj(round(size(proj,1)/2),:,:))
ishow( proj(:,round(size(proj,2)/2),:))
proj_size = size( proj );
cl = class( proj);
fprintf('\nproj: %u %u %u, %s, %f MG, %f GB', proj_size, cl, MB( proj ), GB( proj) )
domain(proj)
fprintf('\n read in %f s = %f min', t, t/60 )



nn = nn + 1;

fn_hdf = sprintf( '%s/hdf/test_%04u.h5', out_path, nn );
fn_bin = sprintf( '%s/binary/test_%04u.raw', out_path, nn );
fn_tif_multi = sprintf( '%s/tif_multi/test_single_image_%04u.tif', out_path, nn );
fn_tif_single = sprintf( '%s/tif_single/run%04u/', out_path, nn );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n\nWRITE TEST' )

%%%%%%%%%%%%%%%%%%%
fprintf( '\n HDF5' )


CheckAndMakePath( fn_tif_single )

tic
h5create( fn_hdf, '/dataset1', proj_size, 'Datatype', cl)
thdfcreate = toc;
fprintf( '\n  create in %f s', thdfcreate )

tic
h5write( fn_hdf, '/dataset1', proj )
thdfwrite = toc;
fprintf( '\n  write in %f s = %f min', thdfwrite, thdfwrite/60 )


%%%%%%%%%%%%%%%%%%%
fprintf( '\n TIFF' )

% fprintf( '\n write multi tiff\n  ' )
% tic
% imwrite( proj(:,:,1), fn_tif_multi )
% for n = 2:proj_size(3)
%     imwrite( proj(:,:,n), fn_tif_multi, 'WriteMode', 'append' )
% end
% %write32bitTIF( fnt, im0 )
% toc
tic
format_single = isa(proj, 'single');
parfor n = 1:proj_size(3)
    fn = sprintf( '%sim_%07u.tif', fn_tif_single, n );
    if format_single
        write32bitTIFfromSingle(fn, proj(:,:,n))
    else
        imwrite( proj(:,:,n), fn , 'tiff');
    end
end
ttifwrite = toc;
fprintf( '\n  write in %f s = %f min', ttifwrite, ttifwrite/60 )

%%%%%%%%%%%%%%%%%%%
fprintf( '\n BINARY' )
tic
    fid = fopen( fn_bin, 'w' );
    fwrite( fid, proj, cl );
    fclose( fid );
tbinwrite = toc;
fprintf( '\n  write in %f s = %f min', tbinwrite, tbinwrite/60 )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n\nREAD TEST' )

%%%%%%%%%%%%%%%%%%%
fprintf( '\n HDF' )
tic
vol = h5read( fn_hdf, '/dataset1' );
thdfread = toc;
fprintf( '\n  read in %f s = %f min', thdfread, thdfread/60 )

%%%%%%%%%%%%%%%%%%%
fprintf( '\n TIFF' )
tic
vol = zeros( [imsize numel( fs )], cl );
imf1 = imfinfo(sprintf( '%sim_%07u.tif', fn_tif_single, 1 ) );
parfor n = 1:num_proj
    fn = sprintf( '%sim_%07u.tif', fn_tif_single, n );
    %proj(:,:,n) = imread( fn, 'tif', 'Info', imf1 )
    vol(:,:,n) = imread( fn, 'tif' )
end
ttifread = toc;
fprintf( '\n  read in %f s = %f min', ttifread, ttifread/60 )

%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n BINARY' )
tic
fid = fopen( fn_bin, 'r' );
[vol, cnt] = fread( fid, [cl '=>' cl]);
fclose( fid );
tbinread = toc;
fprintf( '\n  read in %f s = %f min', tbinread, tbinread/60 )
tic
vol = reshape( vol, proj_size );
tbinreshape = toc;
fprintf( '\n  reshape in %f s = %f min', tbinreshape, tbinreshape/60 )

% TABLE
fprintf( '\n              ' )

fprintf( '\n               ' )
fprintf( '%15s', 'hdf', 'hdf create', 'tif', 'bin')
fprintf( '\n writing / s:  ' )
fprintf( '%15g', thdfwrite, (thdfcreate), ttifwrite, tbinwrite )
fprintf( '\n writing / min:' )
fprintf( '%15g', [thdfwrite, (thdfcreate), ttifwrite, tbinwrite] / 60 )

fprintf( '\n\n               ' )
fprintf( '%15s', 'hdf', 'tif', 'bin', 'reshape')
fprintf( '\n reading / s:   ' )
fprintf( '%15g', thdfread, ttifread, tbinread, tbinreshape )
fprintf( '\n reading / min:' )
fprintf( '%15g', [thdfread, ttifread, tbinread, tbinreshape ] / 60)
fprintf( '\n' )


%%
fprintf( '\n Orthogonal read test (sinogram)' )
inpath =  [out_path '/hdf/'];
fs = dir( [inpath '*.h5'] );
fn = [fs(1).folder filesep fs(1).name];

%% Cube
%fn = '/asap3/petra3/gpfs/p07/2020/data/11010172/scratch_cc/test_io/hdf/cube4000.h5';
vol_size = 2000 + [1 2 3];
fn = [ out_path '/hdf/cube2000.h5'];
fprintf( '\nWrite cube to hdf' )
cl = 'uint16';
tic
h5create( fn, '/dataset1', vol_size, 'Datatype', cl)
h5write( fn, '/dataset1', ones( vol_size, cl) )
thdfwrite = toc;
fprintf( '\n write in %f s = %f min', thdfwrite, thdfwrite/60 )
%%
h5disp( fn )
h5i = h5info( fn );
vol_size = h5i.Datasets.Dataspace.Size;
v = vol_size(1);
h = vol_size(2);
a = vol_size(3);

%% sino
fprintf('Read sino from hdf')
h_start = [round( v / 2 ) 1 1];
h_count = [1 h a];
h_stride = [1 1 1];
fprintf( '\n start : %u %u %u', h_start )
fprintf( '\n count : %u %u %u', h_count )
fprintf( '\n stride: %u %u %u', h_stride )
fprintf( '\n start reading subset' )
tic
hsino = h5read( fn, '/dataset1', h_start, h_count, h_stride );
tortho = toc;
fprintf( '\n read sino in %f s = %f min', tortho, tortho/60 )
%% proj
fprintf('Read proj form hdf')
h_start = [1 1 round( a / 2)];
h_count = [v h 1];
h_stride = [1 1 1];
fprintf( '\n start : %u %u %u', h_start )
fprintf( '\n count : %u %u %u', h_count )
fprintf( '\n stride: %u %u %u', h_stride )
fprintf( '\n start reading subset' )
tic
hproj = h5read( fn, '/dataset1', h_start, h_count, h_stride );
tproj = toc;
fprintf( '\n read proj in %f s = %f min\n', tproj, tproj/60 )
%% 2n ortho
fprintf('Read remaining dimension from hdf')
h_start = [1 round( h / 2) 1];
h_count = [v 1 a];
h_stride = [1 1 1];
fprintf( '\n start : %u %u %u', h_start )
fprintf( '\n count : %u %u %u', h_count )
fprintf( '\n stride: %u %u %u', h_stride )
fprintf( '\n start reading subset' )
tic
hortho = h5read( fn, '/dataset1', h_start, h_count, h_stride );
tortho = toc;
fprintf( '\n read ortho in %f s = %f min\n', tortho, tortho/60 )

fprintf( '\n' )
