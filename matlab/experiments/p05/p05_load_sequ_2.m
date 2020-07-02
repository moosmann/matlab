function [vol, roi] = p05_load_sequ_2( p )
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
steps = assign_from_struct( p , 'steps', [] ); % # of tomos to process, use low number for testing
out_thresh = assign_from_struct( p , 'out_thresh', 0.01 ); % percentage of outliers to be thresholded
roi = assign_from_struct( p , 'roi', [] ); % crop range

% TO DO
% get outlier thresholds from several slices

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
ca
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

if isscalar( steps )
    steps = 1:steps;
end
num_steps = numel( steps );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Read in load sequence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = toc;
fprintf( '\nRead in %u tomograms:', numel( steps ) )
for nn = 1:numel( steps )
    fprintf( '\n %2u : %s.', nn, struct_scans(nn).name)
    p = [proc_path filesep struct_scans(nn).name reco_sub];
    struct_slices = dir( [p filesep '*.tif']);
    num_slices = numel( struct_slices );
    parfor mm = 1:num_slices
        impath = [p filesep struct_slices(mm).name];
        im = imread( impath );
        vol(:,:,mm,nn) = conv8bit( im );
    end
    fprintf(  ' Slices : %u.', num_slices )
    pause(0.01)
end
fprintf( '\n Read data in %.1f s', toc - t);

%% Crop ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty( roi )
    fprintf( '\n Mean along 4th dim' )
    vol_mean4 = mean( vol, 4 );
    
    % Find image horizontal crop range
    fprintf( '\n Mean along 3rd dim' )
    vol_mean4_mean3 = mean( vol_mean4, 3 );
    f = figure('Name', 'INTERACTIVE CROP TOOL: horizontal', 'WindowState', 'maximized' );
    [imc, rect] = imcrop( vol_mean4_mean3, [] );
    close( f );
    rect = round( rect );
    x0 = rect(2);
    x1 = x0 + rect(4);
    y0 = rect(1);
    y1 = y0 + rect(3);
    figure( 'Name', 'cropped image horizontally' )
    imsc( imc )
    axis tight equal
    
    % Find image vertical range
    fprintf( '\n Mean along 2nd dim' )
    vol_mean4_mean2 = squeeze( mean( vol_mean4, 2 ) );
    f = figure('Name', 'INTERACTIVE CROP TOOL: vertical', 'WindowState', 'maximized' );
    [imc, rect] = imcrop( vol_mean4_mean2, [] );
    close( f );
    rect = round( rect );
    z0 = rect(1);
    z1 = z0 + rect(3);
    figure( 'Name', 'cropped image vertically' )
    imsc( imc )
    axis tight equal
    
    roi = [x0 x1 y0 y1 z0 z1];
    fprintf( '\n p.roi = [%u %u %u %u %u %u];', roi )
else
    x0 = roi(1);
    x1 = roi(2);
    y0 = roi(3);
    y1 = roi(4);
    z0 = roi(5);
    z1 = roi(6);
end
% Crop
vol = vol( x0:x1, y0:y1, z0:z1, :);

%% Save 3D binary  tif %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binary_thresh = 145;
fprintf( '\n Save binary 3D multi tif:' )
t = toc;
tif_path = [scan_path filesep 'tif_binary'];
CheckAndMakePath( tif_path )
volb = vol > binary_thresh;
for mm = 1:size( volb, 4 )
    fprintf( ' %u', mm )
    filename = sprintf( '%s/%s_loadSequ_%02u.tif', tif_path, scan_name, mm );
    for nn = 1:size( volb, 3 )
        im = squeeze( volb(:,:,nn,mm) );
        if nn == 1
            imwrite( im, filename, 'Compression','none');
        else
            imwrite( im, filename, 'WriteMode', 'append',  'Compression','none');
        end
    end
end
fprintf( ' Done in %.1f s', toc - t);

%% Save 3D tif %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Save 3D multi tif:' )
t = toc;
tif_path = [scan_path filesep 'tif'];
CheckAndMakePath( tif_path )
for mm = 1:size( vol, 4 )
    fprintf( ' %u', mm )
    filename = sprintf( '%s/%s_loadSequ_%02u.tif', tif_path, scan_name, mm );
    for nn = 1:size( vol, 3 )
        im = squeeze( vol(:,:,nn,mm) );
         if nn == 1
            imwrite( im, filename, 'Compression','none');
        else
            imwrite( im, filename, 'WriteMode', 'append',  'Compression','none');
        end
    end
end
fprintf( ' Done in %.1f s', toc - t);

%% Save animated gif slices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Save animated gifs' )
t = toc;
[sx, sy, sz, ~] = size( vol );
xx = round( sx / 2 );
yy = round( sy / 2 );
zz = round( sz / 2 );
gif_path = [scan_path filesep 'gif'];
CheckAndMakePath( gif_path )
% x
filename = sprintf( '%s/%s_loadSequ_x%04u.gif', gif_path, scan_name, xx );
fprintf( '\n output file: %s', filename)
map = colormap(gray);
for nn = 1:size(vol,4)
    %[A, map] = gray2ind(rot90(squeeze( vol(xx,:,:,nn))));
    A = gray2ind(rot90(squeeze( vol(xx,:,:,nn))));
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
% y
filename = sprintf( '%s/%s_loadSequ_y%04u.gif', gif_path, scan_name, yy );
fprintf( '\n output file: %s', filename)
for nn = 1:size(vol,4)
    A = gray2ind(rot90(squeeze( vol(:,yy,:,nn))));
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
% z
filename = sprintf( '%s/%s_loadSequ_z%04u.gif', gif_path, scan_name, zz );
fprintf( '\n output file: %s', filename)
for nn = 1:size(vol,4)
    A = gray2ind((squeeze( vol(:,:,zz,nn))));
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
fprintf( '\n Done in %.1f s', toc - t);

%% Show %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nimplay( squeeze( vol(xx,:,:,:) ) )
%nimplay( squeeze( vol(:,yy,:,:) ) )
nimplay( cat(2, flipud(permute( squeeze(vol(xx,:,:,:)), [ 2 1 3])), flipud(permute( squeeze(vol(:,yy,:,:)), [ 2 1 3]))), 1, [], 'x-/y-cut sequence' )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Total time elapsed: %.1f s', toc);
fprintf('\nFINISHED\n')