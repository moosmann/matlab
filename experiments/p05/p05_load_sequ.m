function [vol, vol_reg] = p05_load_sequ( proc_path, scan_name, reco_sub, regdir, steps, out_thresh, register, auto_roi_center )
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
if nargin < 1
    proc_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/processed';
end
if nargin < 2
    scan_name = 'syn007_94L_Mg10Gd_8w';
end
if nargin < 3
    reco_sub = '/reco/float_rawBin4';
end
if nargin < 4
    regdir = 'x'; % x or y cut, y sometimes works better if there are more structures to correlate
end
if nargin < 5
    steps = []; % # of tomos to process, use low number for testing
end
if nargin < 6
    out_thresh = 0.01; % percentage of outliers to be thresholded
end
if nargin < 7
    register = 0;
end
if nargin < 8
    auto_roi_center = 0;
end
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
for nn = 1:numel( steps )
    fprintf( '\n %2u : %s.', nn, struct_scans(nn).name)
    p = [proc_path filesep struct_scans(nn).name reco_sub];
    struct_slices = dir( [p filesep '*.tif']);
    num_slices = numel( struct_slices );
    fprintf('\n %u', nn)
    %parfor mm = 1:num_slices
    for mm = 1:num_slices
        impath = [p filesep struct_slices(mm).name];
        vol(:,:,mm,nn) = conv8bit( imread( impath ) );
    end
    fprintf(  ' Slices : %u.', num_slices )
    pause(0.01)
end
fprintf( '\n Read data in %.1f s', toc - t);

%% Auto ROI
if auto_roi_center
    roi_cen = zeros( [3, 2] );
    for nn = 1:3
        im = squeeze( mean( squeeze( vol(:,:,:,1) ), nn ) );
        for mm = 1:2
            im_std = squeeze( std( squeeze( im ), 0, mm ) ) ;
            im_std = SubtractMean( im_std );
            im_std( im_std <= 0 ) = 0;
            roi_cen(nn,mm) = round( CenterOfMass( im_std ) );
        end
    end
    xx = roi_cen(3,1);
    yy = roi_cen(3,2);
    zz = round(( roi_cen(1,2) + roi_cen(2,2) ) / 2);
else
    xx = round( size( vol, 1) / 2);
    yy = round( size( vol, 2) / 2);
    zz = round( size( vol, 3) / 2);
end
fprintf( '\nroi center before registering: x,y,z = %u,%u,%u',xx,yy,zz)

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
        xx = roi_cen(3,1);
        yy = roi_cen(3,2);
        zz = round( ( roi_cen(1,2) + roi_cen(2,2) ) / 2 );
    else
        xx = round( size( vol_reg, 1) / 2);
        yy = round( size( vol_reg, 2) / 2);
        zz = round( size( vol_reg, 3) / 2);
    end
    fprintf( '\n roi center after registering: x,y,z = %u,%u,%u',xx,yy,zz)
    
    
else
    fprintf( '\n No registerion of slice ' )
    vol_reg = vol;
end

%% Save
%nimplay(cat(3,im1(s:end,:),im2(1:end-s+1,:)))
%nimplay( squeeze(vol_reg(xx,:,:,:)));

% Slice: Animated gif
t = toc;
fprintf( '\n Save animated gifs' )

CheckAndMakePath( scan_path )
filename = sprintf( '%s/%s_loadSequ_x%04u.gif', scan_path, scan_name, xx );
fprintf( '\n output file: %s', filename)

map = colormap(gray);
for nn = 1:size(vol_reg,4)
    %[A, map] = gray2ind(rot90(squeeze( vol_reg(xx,:,:,nn))));
    A = gray2ind(rot90(squeeze( vol_reg(xx,:,:,nn))));
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

filename = sprintf( '%s/%s_loadSequ_y%04u.gif', scan_path, scan_name, yy );
fprintf( '\n output file: %s', filename)

for nn = 1:size(vol_reg,4)
    A = gray2ind(rot90(squeeze( vol_reg(:,yy,:,nn))));
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

filename = sprintf( '%s/%s_loadSequ_z%04u.gif', scan_path, scan_name, zz );
fprintf( '\n output file: %s', filename)

for nn = 1:size(vol_reg,4)
    A = gray2ind((squeeze( vol_reg(:,:,zz,nn))));
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
fprintf( '\n Done in %.1f s', toc - t);

% Save 4D matlab array
% t = toc;
% fprintf( '\n Save 4D matlab volume.' )
% filename = sprintf( '%s/%s_loadSequ_4D.mat', scan_path, scan_name );
% fprintf( '\n output file: %s', filename)
% save( filename, 'vol_reg' , '-v7.3', '-nocompression')
% fprintf( '\n Done in %.1f s', toc - t);

%% Show %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nimplay( squeeze( vol_reg(xx,:,:,:) ) )
%nimplay( squeeze( vol_reg(:,yy,:,:) ) )
nimplay( cat(2, flipud(permute( squeeze(vol_reg(xx,:,:,:)), [ 2 1 3])), flipud(permute( squeeze(vol_reg(:,yy,:,:)), [ 2 1 3]))), 1, [], 'x-/y-cut sequence' )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Total time elapsed: %.1f s', toc);
fprintf('\nFINISHED\n')