function [vol, vol_reg] = p05_load_sequ( p )
% Create images sequences from 4D load tomography data.
%
% Script to process in 4D image squences e.g. from load tomography. Read in
% tomo data sets, convert to 8bit scaled by given thresholds, and save load
% sequences cut centrally along x and y as animated gifs. After processing
% registered 4D volume is available as 4D array in the Matlab workspace.
% (4D array is not saved since saving the Matlab array is more time
% consuming than processing the sequence.)
%
% proc_name : string, name of sequence without trailing indices and
%   underscore
% reco_sub : string, subfolder to 'reco' folder
% out_thresh : scalar, < 1, percentage of outliers to be filtered before
%   8bit conversion

%% DEFAULT ARGUMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proc_path = assign_from_struct( p , 'proc_path', '/asap3/petra3/gpfs/p05/2017/data/11004016/processed' );
scan_name = assign_from_struct( p, 'scan_name', 'syn007_94L_Mg10Gd_8w' );
reco_sub = assign_from_struct( p , 'reco_sub', 'float_rawBin2' );
regdir = assign_from_struct( p , 'regdir', 'x' ); % x or y cut, y sometimes works better if there are more structures to correlate
steps = assign_from_struct( p , 'steps', [] ); % # of tomos to process, use low number for testing
out_thresh = assign_from_struct( p , 'out_thresh', 0.01 ); % percentage of outliers to be thresholded
register = assign_from_struct( p , 'register', 0 );
auto_roi_center = assign_from_struct( p , 'auto_roi_center', 0 );
crop_roi = assign_from_struct( p , 'crop_roi', [] );
barcol = assign_from_struct( p , 'barcol', 'white' );
voxel_size = assign_from_struct( p , 'voxel_size', [] );
proj_type = assign_from_struct( p , 'proj_type', 'max' );
% TO DO
% get outlier thresholds from several slices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

fprintf( '\n LOAD SEQUENCE PROCESSING')

scan_path = [proc_path filesep scan_name];
CheckAndMakePath( scan_path )

fprintf( '\n scan path : %s', scan_path)
fprintf( '\n scan name : %s', scan_name)

% scans to process
struct_scans = dir([scan_path '_*']);
if numel( struct_scans ) == 0
    error( 'No recos found matching folder pattern: \n%s\n', scan_path )
end

if isempty( steps )
    steps = 1:numel( struct_scans);
end

if ~strcmp( reco_sub(1), '/' )
    reco_sub = [ '/' reco_sub];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isscalar( steps )
    steps = 1:steps;
end
num_steps = numel( steps );

% Preallocation
fprintf( '\n Allocation of 4D volume. Size:' )
p = [proc_path filesep struct_scans(3).name reco_sub];
struct_slices = dir( [p filesep '*.tif']);
num_slices = numel( struct_slices );
tinfo = imfinfo( [struct_slices(1).folder filesep struct_slices( floor( num_slices / 2 )).name]);
im = imread( [struct_slices(1).folder filesep struct_slices( floor( num_slices / 2 )).name]);
vol = zeros( [tinfo.Width tinfo.Height num_slices num_steps ] , 'uint8');
fprintf( ' %u ', size( vol) )
fprintf( '. Memory: %f GiB ', GB( vol ))

% Outlier threshold
im_sorted = sort( im(:) );
ind = round( out_thresh * numel( im_sorted ) );
im_low = im_sorted(ind);
im_high = im_sorted(end-ind);
fprintf( '\n image : min, max = %f, %f', im_sorted(1), im_sorted(end) )
fprintf( '\n %g%% thresholds: low, high = %f, %f', out_thresh, im_low, im_high )

conv8bit = @(im) uint8( (2^8 - 1) * ( im - im_low) / ( im_high - im_low ) );
imc = conv8bit( im );

% Figure to check conversion scaling
figure( 'Name', 'Check thresholding for 8bit conversion')
subplot( 1, 2, 1 )
imsc(im);
title(sprintf('full range:32bit single float'))
axis equal tight
subplot( 1, 2, 2 )
imsc(imc);
title(sprintf('outlier thresholding: 8bit uint'))
axis equal tight
drawnow

%% Read in load sequence
t = toc;
fprintf( '\nRead in %u tomograms:', numel( steps ) )
for ss = 1:numel( steps )
    nn = steps(ss);
    fprintf( '\n %2u : %s.', nn, struct_scans(nn).name)
    p = [proc_path filesep struct_scans(nn).name reco_sub];
    struct_slices = dir( [p filesep '*.tif']);
    num_slices = numel( struct_slices );
    parfor mm = 1:num_slices        
        impath = [p filesep struct_slices(mm).name];
        im = conv8bit( imread( impath ) ) ;
        im = rot90(im);
        vol(:,:,mm,ss) = im;
    end
    fprintf(  ' Slices : %u.', num_slices )
    pause(0.01)
end
[sx, sy, sz, ~] = size( vol );
fprintf( '\n Read data in %.1f s', toc - t);

switch proj_type
    case 'mean'
        vol_proj = mean( vol,4 );
    case 'median'
        vol_proj = median( vol,4 );
    case 'max'
        vol_proj = max( vol, [],4 );
    case 'min'
        vol_proj = min( vol, [],4 );
end

%% Vertical ROI
im1max = normat(  double( squeeze(max( vol_proj,[],1)) ) );
im1min = normat( -double( squeeze(min( vol_proj,[],1)) ) );
im2max = normat(  double( squeeze(max( vol_proj,[],2)) ) );
im2min = normat( -double( squeeze(min( vol_proj,[],2)) ) );
figure('Name', 'Extrema projections' )
subplot( 2, 2, 1 )
imsc( im1max )
axis equal tight
subplot( 2, 2, 2 )
imsc( im1min )
axis equal tight
subplot( 2, 2, 3 )
imsc( im2max )
axis equal tight
subplot( 2, 2, 4 )
imsc( im2min )
axis equal tight
l = mean( im1max + im1min + im2max + im2min, 1 );
lm = mean( l );
for zz = 1:sz
    z0 = zz;
    if l(z0) > lm
        break;
    end
end
for zz = sz:-1:z0+1
    z1 = zz;
    if l(z1) > lm
        break;
    end
end
dz = 0.1* (z1-z0);
z0 = max( 1, z0 - dz );
z1 = min( sz, z1 + dz );
fprintf( '\n vertical ROI: %u %u', z0, z1 )

%% Auto ROI
if auto_roi_center
    roi_cen = zeros( [3, 2] );
    for nn = 1:3
        im = squeeze( mean( squeeze( vol_proj ), nn ) );
        for mm = 1:2
            im_std = squeeze( std( squeeze( im ), 0, mm ) ) ;
            im_std = SubtractMean( im_std );
            im_std( im_std <= 0 ) = 0;
            roi_cen(nn,mm) = round( CenterOfMass( im_std ) );
        end
    end
    xx = roi_cen(3,2);
    yy = roi_cen(3,1);
    zz = round(( roi_cen(1,2) + roi_cen(2,2) ) / 2);
else
    xx = round( size( vol, 1) / 2);
    yy = round( size( vol, 2) / 2);
    zz = round( size( vol, 3) / 2);
end
fprintf( '\n roi center : x,y,z = %u,%u,%u',xx,yy,zz)
if ~isempty( crop_roi ) && crop_roi
    cx = round( crop_roi(1) / 2 );
    cy = round( crop_roi(2) / 2 );
    x0 = max( 1, xx - cx );
    x1 = min( sx , xx + cx);
    y0 = max( 1, yy - cy);
    y1 = min( sy , yy + cy);
    vol= vol(x0:x1,y0:y1,z0:z1,:);
    [sx, sy, sz, ~] = size( vol );
    fprintf( '\n ROI shape : %u %u %u', sx, sy, sz );
end
figure( 'Name', 'vol roi' )
imsc( vol(:,:,round(sz/2),1))
axis equal tight

%% Registering
if register
    t = toc;
    fprintf( '\n Registering slices ' )
    shift = zeros( [num_steps, 1]);
    for nn = 1:num_steps - 1
        if strcmp( regdir, 'x')
            im1 = squeeze( vol(xx,100:end-100,100:end-100,nn));
            im2 = squeeze( vol(xx,100:end-100,100:end-100,nn + 1));
        elseif strcmp( regdir, 'y' )
            im1 = squeeze( vol(100:end-100,yy,100:end-100,nn));
            im2 = squeeze( vol(100:end-100,yy,100:end-100,nn + 1));
        elseif strcmp( regdir, 'z' )
            im1 = squeeze( vol(100:end-100,100:end-100,zz,nn));
            im2 = squeeze( vol(100:end-100,100:end-100,zz,nn + 1));
        end
        out = ImageCorrelation(im1,im2);
        shift(nn+1) = round( out.shift2);
    end
    z0 = abs( cumsum( shift ) );
    z1 = max(z0) - z0;
    fprintf( '\n frame-to-frame shifts : ' )
    fprintf( '%g ', shift)
    fprintf( '\n global shifts : ' )
    fprintf( '%g ', z0)
    
    fprintf( '\n Cropping volumes. Original size : %u x %u x %u x %u.', size( vol) )
    vol_reg = zeros( [size(vol,1), size(vol,2), size(vol,3) - max(z0), num_steps], 'uint8' );
    fprintf( '\n New size : %u x %u x %u x %u ', size( vol_reg) )
    for nn = 1:num_steps
        vol_reg(:,:,:,nn) = vol(:,:,1+z0(nn):end-z1(nn),nn);
        %disp( size( vol(:,:,1+z0(nn):end-z1(nn),nn) ) )
        
    end
    fprintf( 'Done in %.1f s', toc - t);
else
    fprintf( '\n No registration of slices ' )
    vol_reg = vol;
end

%% Auto ROI
if auto_roi_center
    roi_cen = zeros( [3, 2] );
    for nn = 1:3
        im = squeeze( mean( squeeze( vol_reg(:,:,:,1) ), nn ) );
        for mm = 1:2
            im_std = squeeze( std( squeeze( im ), 0, mm ) ) ;
            im_std = SubtractMean( im_std );
            im_std( im_std <= 0 ) = 0;
            roi_cen(nn,mm) = round( CenterOfMass( im_std ) );
        end
    end
    xx = roi_cen(3,2);
    yy = roi_cen(3,1);
    zz = round( ( roi_cen(1,2) + roi_cen(2,2) ) / 2 );
else
    xx = round( size( vol_reg, 1) / 2);
    yy = round( size( vol_reg, 2) / 2);
    zz = round( size( vol_reg, 3) / 2);
end
fprintf( '\n roi center: x,y,z = %u,%u,%u',xx,yy,zz)

CheckAndMakePath( scan_path )

%% Slice: Animated gif
gif_path = [ scan_path filesep 'gif' ];
CheckAndMakePath( gif_path )
t = toc;
fprintf( '\n Save animated gifs' )
map = colormap(gray);

% gif x
% Scalebar coordinates
A = gray2ind(rot90(squeeze( vol_reg(xx,:,:,1)), -1 ));
[dx,dy] = size( A );
barWidth = dy / 10;
fac = 1;
while round( barWidth / fac ) > 0
    barWidth = fac * round( barWidth / fac );
    fac = fac * 10;
end
margin = max( 10, round( min( dx, dy ) / 20 ) );
barHeight = max( 2, round( dx / 50 ) );
barx = 1+dx-barHeight-margin:dx-margin;
bary = 1+dy-barWidth-margin:dy-margin;
switch barcol
    case 'white'
        barval = max(A(:));
    case 'black'
        barval = max(A(:));
end
fprintf( '\n scalebar width : %u', barWidth )
fprintf( '\n scalebar heigth : %u', barHeight )
fprintf( '\n scalebar margin : %u', margin )
fprintf( '\n scalebar grayvalue : %g', barval )
barstr = sprintf( 'barH%uW%u', barHeight, barWidth );
if ~isempty( voxel_size )
    barstr = sprintf( '%s_voxel%.2f', barstr, voxel_size*1e6 );
    barstr = regexprep( barstr, '\.', 'p' );
end
fprintf( '\n scalebar string : %s', barstr )
filename = sprintf( '%s/%s_loadSequ_x%04u_%s.gif', gif_path, scan_name, xx, barstr );
fprintf( '\n output file: %s', filename)
for nn = 1:size(vol_reg,4)
    A = gray2ind(rot90(squeeze( vol_reg(xx,:,:,nn))));
    A(barx,bary) = barval;
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

% gif y
% Scalebar coordinates
A = gray2ind(rot90(squeeze( vol_reg(:,yy,:,1)), -1));
[dx,dy] = size( A );
barWidth = dy / 10;
fac = 1;
while round( barWidth / fac ) > 0
    barWidth = fac * round( barWidth / fac );
    fac = fac * 10;
end
margin = max( 10, round( min( dx, dy ) / 20 ) );
barHeight = max( 2, round( dx / 50 ) );
barx = 1+dx-barHeight-margin:dx-margin;
bary = 1+dy-barWidth-margin:dy-margin;
switch barcol
    case 'white'
        barval = max(A(:));
    case 'black'
        barval = max(A(:));
end
fprintf( '\n scalebar width : %u', barWidth )
fprintf( '\n scalebar heigth : %u', barHeight )
fprintf( '\n scalebar margin : %u', margin )
fprintf( '\n scalebar grayvalue : %g', barval )
barstr = sprintf( 'barH%uW%u', barHeight, barWidth );
if ~isempty( voxel_size )
    barstr = sprintf( '%s_voxel%.2f', barstr, voxel_size*1e6 );
    barstr = regexprep( barstr, '\.', 'p' );
end
fprintf( '\n scalebar string : %s', barstr )
filename = sprintf( '%s/%s_loadSequ_y%04u_%s.gif', gif_path, scan_name, yy, barstr );
fprintf( '\n output file: %s', filename)
for nn = 1:size(vol_reg,4)
    A = gray2ind(rot90(squeeze( vol_reg(:,yy,:,nn))));
    A(barx,bary) = barval;
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

% gif z
% Scalebar coordinates
A = gray2ind((squeeze( vol_reg(:,:,zz,1))));
[dx,dy] = size( A );
barWidth = dy / 10;
fac = 1;
while round( barWidth / fac ) > 0
    barWidth = fac * round( barWidth / fac );
    fac = fac * 10;
end
margin = max( 10, round( min( dx, dy ) / 20 ) );
barHeight = max( 2, round( dx / 50 ) );
barx = 1+dx-barHeight-margin:dx-margin;
bary = 1+dy-barWidth-margin:dy-margin;
switch barcol
    case 'white'
        barval = max(A(:));
    case 'black'
        barval = max(A(:));
end
fprintf( '\n scalebar width : %u', barWidth )
fprintf( '\n scalebar heigth : %u', barHeight )
fprintf( '\n scalebar margin : %u', margin )
fprintf( '\n scalebar grayvalue : %g', barval )
barstr = sprintf( 'barH%uW%u', barHeight, barWidth );
if ~isempty( voxel_size )
    barstr = sprintf( '%s_voxel%.2f', barstr, voxel_size*1e6 );
    barstr = regexprep( barstr, '\.', 'p' );
end
fprintf( '\n scalebar string : %s', barstr )
filename = sprintf( '%s/%s_loadSequ_z%04u_%s.gif', gif_path, scan_name, zz, barstr );
fprintf( '\n output file: %s', filename)
for nn = 1:size(vol_reg,4)
    A = gray2ind((squeeze( vol_reg(:,:,zz,nn))));
    A(barx,bary) = barval;
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
fprintf( '\n Done in %.1f s', toc - t);

%% Save 3D tif z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Save 3D multi tifs z .' )
t = toc;
tif_path = [scan_path filesep 'tif'];
CheckAndMakePath( tif_path )
for mm = 1:size( vol_reg, 4 )
    filename = sprintf( '%s/%s_loadSequ_z_%02u.tif', tif_path, scan_name, mm );
    for nn = 1:size( vol_reg, 3 )
        im = squeeze( vol_reg(:,:,nn,mm) );
        if nn == 1
            imwrite( im, filename, 'Compression','none');
        else
            imwrite( im, filename, 'WriteMode', 'append',  'Compression','none');
        end
    end
end
fprintf( ' Done in %.1f s', toc - t);

%% Save 3D tif x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Save 3D multi tifs x slicing.' )
t = toc;
tif_path = [scan_path filesep 'tif'];
CheckAndMakePath( tif_path )
for mm = 1:size( vol_reg, 4 )
    filename = sprintf( '%s/%s_loadSequ_x_%02u.tif', tif_path, scan_name, mm );
    for nn = 1:size( vol_reg, 3 )
        im = squeeze( vol_reg(:,:,nn,mm) );
        if nn == 1
            imwrite( im, filename, 'Compression','none');
        else
            imwrite( im, filename, 'WriteMode', 'append',  'Compression','none');
        end
    end
end
fprintf( ' Done in %.1f s', toc - t);
%% Save ROI tif
if ~isempty( crop_roi )
    for mm = 1:size( vol_reg, 4 )
        filename = sprintf( '%s/%s_loadSequ_%02u.tif', scan_path, scan_name, mm );
        for nn = 1:size( vol_reg, 3 )
            im = squeeze( vol_reg(:,:,nn,mm) );
            if nn == 1
                imwrite( im, filename, 'Compression','none');
            else
                imwrite( im, filename, 'WriteMode', 'append',  'Compression','none');
            end
        end
    end
end

%% Show %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nimplay( squeeze( vol_reg(xx,:,:,:) ) )
%nimplay( squeeze( vol_reg(:,yy,:,:) ) )
nimplay( cat(2, flipud(permute( squeeze(vol_reg(xx,:,:,:)), [ 2 1 3])), flipud(permute( squeeze(vol_reg(:,yy,:,:)), [ 2 1 3]))), 1, [], 'x-/y-cut sequence' )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Total time elapsed: %.1f s', toc);
fprintf('\nFINISHED\n')