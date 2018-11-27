%vol_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn151_58L_Mg_12_001/reco/uint8_rawBin2_recoBin2';
vol_path = '/asap3/petra3/gpfs/p05/2018/data/11004263/processed/syn011_90R_Mg5Gd_4w_a/reco/float_rawBin2_recoBin2';

% Read volume
fprintf( '\n Read data\n')
StepSize_or_VecOfImagesToRead = 1;
FilenamePattern = '*.tif';
raw_im_shape = []; % for raw files only
if ~exist( 'vol', 'var' )
    vol = read_images_to_stack( vol_path, StepSize_or_VecOfImagesToRead, FilenamePattern, raw_im_shape );
end

% sum up images along all 3 dimensions
vol_mean1 = squeeze( mean( vol, 1 ) );
vol_mean2 = squeeze( mean( vol, 2 ) );
vol_mean3 = squeeze( mean( vol, 3 ) );

%% Histogram..
nbins = 100;
h = histogram( vol(:), nbins );

%% Show summed images
if exist( 'f1', 'var' ) && isvalid( f1 )
    figure( f1 )
else
    f1 = figure( 'Name', 'Projections' );
end

subplot(1,3,1)
imsc( vol_mean1)
title( 'mean 1')
axis tight equal
colorbar

subplot(1,3,2)
imsc( vol_mean2)
title( 'mean 2')
axis tight equal
colorbar

subplot(1,3,3)
imsc( vol_mean3 )
title( 'mean 3')
axis tight equal
colorbar

%% Sum projections
vol_mean3_mean1 = mean( vol_mean3, 1 );
vol_mean3_mean2 = mean( vol_mean3, 2 );

%% Indices where to crop
threshold = 0; % Check. maybe it hast to be set by 0.001*max(vol_mean3_mean1) or similar
m = vol_mean3_mean1 > threshold;
% coordinate vector
y = 1:size( vol, 2 ); % CHECK
n = y(m);
% min/max index
y1 = n(1);
y2 = n(end);

m = vol_mean3_mean2 > threshold;
% coordinate vector
x = 1:size( vol, 1 ); % CHECK
n = x(m);
% min/max index
x1 = n(1);
x2 = n(end);

%% Show projections and indeces
if exist( 'f2', 'var' ) && isvalid( f2 )
    figure( f2 )
else
    f2 = figure( 'Name', 'Summed projections' );
end

subplot( 2, 3, 1 )
plot( vol_mean3_mean1 )
title( 'y: mean 3 mean 1' )
hold on
%plot( y1, 0, 'r*', y2, 0, 'r*' )
plot( y1, threshold, 'MarkerIndices', 10 )
plot( y2, threshold, 'MarkerIndices', 10 )
axis tight

subplot( 2, 3, 4 )
plot( vol_mean3_mean2 )
hold on
plot( x1, 0, 'Marker','o')
title( 'x: mean 3 mean 2' )
axis tight

%% Show cropped slices
if exist( 'f3', 'var' ) && isvalid( f3 )
    figure( f3 )
else
    f3 = figure( 'Name', 'Cropped projections' );
end

subplot( 1, 3, 1 )
imsc( vol_mean3( x1:x2, y1:y2 ) );
title( 'xy (mean3) cropped ' )
axis tight equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Finished\n')