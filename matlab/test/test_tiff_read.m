%% read tiff
clear all
fprintf( '\nTIFF READ TEST\n' )
warning( 'off', 'imageio:tifftagsread:expectedTagDataFormat'  )

raw_path = '/asap3/petra3/gpfs/p07/2020/data/11010172/raw/swerim_21_12_oh_a/tiff0000/*tif';
fs = dir( raw_path );
fn = @(n) [fs(n).folder filesep fs(n).name];


tic
imf = imfinfo( fn(1) );
toc;

tic
imf = imfinfo( fn(2) );
toc;

tic
imf = imfinfo( fn(3), 'tif' );
toc

tic
imf = imfinfo( fn(4), 'tif' );
toc

tic
imf = imfinfo( fn(5) );
toc;

tic
imf = imfinfo( fn(6), 'tif' );
toc



tic;
im3 = imread( fn(7), 'tif' );
toc


tic;
im4 = imread( fn(8), 'tif', 'Info', imf );
toc

tic;
im5 = imread( fn(9), 'tif' );
toc


tic;
im6 = imread( fn(10), 'tif', 'Info', imf );
toc

