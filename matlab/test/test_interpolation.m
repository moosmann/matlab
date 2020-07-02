%out_path = '/asap3/petra3/gpfs/p07/2020/data/11010172/scratch_cc/test_interpolation';
%CheckAndMakePath( out_path );

raw_path = '/asap3/petra3/gpfs/p07/2020/data/11010172/raw/swerim_21_12_oh_a/tiff0000/';

%% read tiff
fprintf( '\nINTERPOLATIN TEST' )

fn = '/asap3/petra3/gpfs/p07/2020/data/11010172/raw/swerim_21_12_oh_a/tiff0000/swerim_21_12_oh_a_0001045_img.tif';
warning( 'off', 'imageio:tifftagsread:expectedTagDataFormat'  )

im0 = imread( fn, 'tif' );
im1 = single( im0' );

%%
tgpu = 0;
r = round( 13 * rand( 1) );
x  = single(2:7920 - r);
im = im1(x,:);

xq = single( x(1) + 0.33:1:x(end) ); 
fprintf( '\nsize : %u %u', size( im ) )


fprintf( '\nGPU create arrays' )
xg = gpuArray(x);
xqg = gpuArray( xq );
tic;
img = gpuArray( im );
t = toc;
fprintf( ' %f s', t)
tgpu = tgpu + t;

fprintf( '\nGPU interpolation' )
tic;
vqg = interp1( xg, img, xqg, 'linear' );
t = toc;
fprintf( ' %f s', t)
tgpu = tgpu + t;

tic;
fprintf( '\nGPU gather       ' )
vqgg = gather( vq );
t = toc;
fprintf( ' %f s', t)
tgpu = tgpu + t;

fprintf( '\nGPU total :       %f s', tgpu )


fprintf( '\nCPU interpolation' )
tic;
vq = interp1( x, im, xq, 'linear' );
tcpu = toc;
fprintf( ' %f s', tcpu)

fprintf( '\n (CPU time) / (GPU total time): %f', tcpu / tgpu )
% 
% 
% tic;
% fprintf( '\nSorting GPU array' )
% imgs = sort( img(:) );
% t = toc;
% fprintf( ' %f s', t)
fprintf( '\n' )