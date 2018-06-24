% Create images sequences from 4D load tomography data.
%
% Script to process in 4D image squences e.g. from load tomography. Read in
% tomo data sets, convert to 8bit scaled by given thresholds, and save load
% sequences cut centrally along x and y as animated gifs. After processing
% registered 4D volume is available as 4D array in the Matlab workspace.
% (4D array is not saved since saving the Matlab array is more time
% consuming than processing the sequence.) 
%
% scan_name : string, name of sequence without trailing indices and
%   underscore
% reco_sub : string, subfolder to 'reco' folder
% out_thresh : scalar, < 1, percentage of outliers to be filtered before
%   8bit conversion

tic
close all hidden

% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn_max = []; % # of tomos to process, use low number for testing
regdir = 'x'; % x or y cut, y sometimes works better if there are more structures to correlate

proc_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/processed';
%scan_name = 'syn002_6L_PEEK_4w'; 
%scan_name = 'syn003_92L_Mg10Gd_4w';
%scan_name = 'syn004_84L_Mg10Gd_4w'; regdir = 'y';
%scan_name = 'syn005_81L_Mg5Gd_8w'; 
%scan_name = 'syn006_75R_Mg10Gd_8w';
scan_name = 'syn008_76R_Mg10Gd_8w';
%scan_name = 'syn009_32R_PEEK_8w';
%scan_name = 'syn010_19R_PEEK_4w';
%scan_name = 'syn011_14R_PEEK_4w'; regdir = 'y';
%scan_name = 'syn012_79L_Mg5Gd_8w'; % empty scans [1:3 14:20];

reco_sub = 'float_rawBin4';
out_thresh = 0.01; % percentage of outliers to be thresholded 

% TO DO
% get outlier thresholds from severl slices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n LOAD SEQUENCE PROCESSING')

scan_path = [proc_path filesep scan_name];
CheckAndMakePath( scan_path )

fprintf( '\n scan path : %s', scan_path)
fprintf( '\n scan name : %s', scan_name)

% scans to process
struct_scans = dir([scan_path '_*']);

if isempty( nn_max )
    nn_max = numel( struct_scans);
end

% Preallocation
fprintf( '\n Allocation of 4D volume. Size:' )
p = [proc_path filesep struct_scans(3).name '/reco/' reco_sub];
struct_slices = dir( [p filesep '*.tif']);
num_slices = numel( struct_slices );
tinfo = imfinfo( [struct_slices(1).folder filesep struct_slices( floor( num_slices / 2 )).name]);
im = imread( [struct_slices(1).folder filesep struct_slices( floor( num_slices / 2 )).name]);
vol = zeros( [tinfo.Width tinfo.Height num_slices nn_max] , 'uint8');
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
title(sprintf('full dynamic range'))
axis equal tight
subplot( 1, 2, 2 )
imsc(imc);
title(sprintf('after outlier thresholding'))
axis equal tight
drawnow

%% Read in load sequence
t = toc;
for nn = 1:nn_max
    fprintf( '\n %2u : %s.', nn, struct_scans(nn).name)
    p = [proc_path filesep struct_scans(nn).name '/reco/' reco_sub];
    struct_slices = dir( [p filesep '*.tif']);
    num_slices = numel( struct_slices );
    parfor mm = 1:num_slices
        impath = [p filesep struct_slices(mm).name];
        vol(:,:,mm,nn) = conv8bit( imread( impath ) );
    end
    fprintf(  ' Slices : %u.', num_slices )
    pause(0.01)
end
fprintf( '\n Read data in %.1f s', toc - t);

%% Registering
t = toc;
fprintf( '\n Registering slices ' )
shift = zeros(nn_max,1);
for nn = 1:nn_max-1
    if strcmp( regdir, 'x')
        im1 = squeeze( vol(round(size(vol,1)/2),100:end-100,100:end-100,nn));
        im2 = squeeze( vol(round(size(vol,1)/2),100:end-100,100:end-100,nn + 1));
    elseif strcmp( regdir, 'y' )
        im1 = squeeze( vol(100:end-100,round(size(vol,2)/2),100:end-100,nn));
        im2 = squeeze( vol(100:end-100,round(size(vol,2)/2),100:end-100,nn + 1));
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

fprintf( '\n Cropping volumes. ' )
vol_reg = zeros( [size(vol,1), size(vol,2), size(vol,3) - max(z0), nn_max], 'uint8' );
for nn = 1:nn_max
    vol_reg(:,:,:,nn) = vol(:,:,1+z0(nn):end-z1(nn),nn);
    %disp( size( vol(:,:,1+z0(nn):end-z1(nn),nn) ) )
    
end
fprintf( 'Done in %.1f s', toc - t);

%% Save
%nimplay(cat(3,im1(s:end,:),im2(1:end-s+1,:)))
%nimplay( squeeze(vol_reg(xx,:,:,:)));

% Slice: Animated gif
t = toc;
fprintf( '\n Save animated gifs' )

xx = round( size( vol_reg, 1) / 2);

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

yy = round( size( vol_reg, 2) / 2);

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

zz = round( size( vol_reg, 3) / 2);

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
nimplay( cat(2, flipud(permute( squeeze(vol_reg(xx,:,:,:)), [ 2 1 3])), flipud(permute( squeeze(vol_reg(:,yy,:,:)), [ 2 1 3]))))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Total time elapsed: %.1f s', toc);
fprintf('\n')
