
flatcor_path = '/asap3/petra3/gpfs/p05/2020/data/11009667/processed/embl_056_200924_dist_1400_zshift_0p2_010/flat_corrected/rawBin2/';
ds_att = dir( [flatcor_path filesep '*.tif'] );

% Intensity / attenuation
filename = sprintf( '%s/%s', flatcor_path,  ds_att(round( numel( ds_att) / 2 )).name );
int = imread( filename );
int = SubtractMean( int );
[N01, N02] = size( int );

if ~exist( 'RECT', 'var')
    [N1c,N2c,int_roi,RECT] = imcrop(int, []);
end

% FFT
int_fft = fft2( padarray( int, size(int), 'both', 'post' ) );
int_roi_fft = fft2( padarray( int_roi, size(int), 'both', 'post' ) );

% Resample FFT ROI on FFT full
[n1, n2] = size( int_roi_fft );
[N1, N2] = size( int_fft );
X = 1:n2;
Y = 1:n1;
V = int_roi_fft;
Xq = 1:(n2 - 1)/(N2-1):n2;
Yq = 1:(n1 - 1)/(N1-1):n1;
int_roi_fft_res = interp2( X, Y', V, Xq, Yq' );

intf =  real( ifft2( int_fft - int_roi_fft_res, 'symmetric' ) );
intf = intf(1:N01,1:N02);


if exist( 'h1' , 'var' ) && isvalid( h1 )
    figure(h1)
else
    h1 = figure( 'Name', 'ROI, fft(ROI)' );
end

subplot(2,1,1)
imsc( int_roi )
title(sprintf('int roi'))
colorbar
axis equal tight

subplot(2,1,2)
im =  normat( 10^-5 + abs( fftshift( fft2( padarray( SubtractMean(int_roi), [size(int)], 'both', 'post') ) ) ) ) ;
imsc( im )
title(sprintf('int roi fft'))
colorbar
axis equal tight


if exist( 'h2' , 'var' ) && isvalid( h2 )
    figure(h2)
else
    h2 = figure( 'Name', 'Filter image' );
end

subplot(2,1,1)
imsc( int )
title(sprintf('int'))
colorbar
axis equal tight

subplot(2,1,2)
imsc( intf )
title(sprintf('int filtered'))
colorbar
axis equal tight

if exist( 'h3' , 'var' ) && isvalid( h3 )
    figure(h3)
else
    h3 = figure( 'Name', 'Filter image' );
end
subplot(1,1,1)
imsc( abs(int - intf ))
title(sprintf('difference'))
colorbar
axis equal tight

