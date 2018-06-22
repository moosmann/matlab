% Script to process in 4D image squences e.g. from load tomography. Read in
% tomo data sets as 8bit scaled to given thresholds.
%clear all
tic
close all hidden

% PARAMETERS
nn_max = []; % # of tomos to process
proc_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/processed';
regdir = 'x';

%scan_name = 'syn002_6L_PEEK_4w'; 
%scan_name = 'syn003_92L_Mg10Gd_4w';
%scan_name = 'syn004_84L_Mg10Gd_4w'; regdir = 'y';
%scan_name = 'syn005_81L_Mg5Gd_8w'; 
scan_name = 'syn006_75R_Mg10Gd_8w';

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
v = zeros( [tinfo.Width tinfo.Height num_slices nn_max] , 'uint8');
fprintf( ' %u ', size( v) )
fprintf( '. Memory: %f GiB ', GB( v ))

% Outlier threshold
im_sorted = sort( im(:) );
ind = round( out_thresh * numel( im_sorted ) );
im_low = im_sorted(ind);
im_high = im_sorted(end-ind);
fprintf( '\n image : min, max = %f, %f', im_sorted(1), im_sorted(end) )
fprintf( '\n %g%% thresholds: low, high = %f, %f', out_thresh, im_low, im_high )

conv8bit = @(im) uint8( (2^8 - 1) * ( im - im_low) / ( im_high - im_low ) );
imc = conv8bit( im );

figure( 'Name', 'Conversion scaling')
subplot( 1, 2, 1 )
imsc(im);
title(sprintf('full range image'))        
axis equal tight
subplot( 1, 2, 2 )
imsc(imc);
title(sprintf('scaled image'))        
axis equal tight

%% Read in load sequence
t = toc;
for nn = 1:nn_max
    fprintf( '\n %2u : %s.', nn, struct_scans(nn).name)
    p = [proc_path filesep struct_scans(nn).name '/reco/' reco_sub];
    struct_slices = dir( [p filesep '*.tif']);
    num_slices = numel( struct_slices );
    parfor mm = 1:num_slices
        impath = [p filesep struct_slices(mm).name];
        v(:,:,mm,nn) = conv8bit( imread( impath ) );
    end
    fprintf(  ' Slices : %u.', num_slices )
end
fprintf( '\n Read data in %.1f s', toc - t);

%% Registering
t = toc;
fprintf( '\n Registering slices ' )
shift = zeros(nn_max,1);
for nn = 1:nn_max-1
    if strcmp( regdir, 'x')
        im1 = squeeze( v(round(size(v,1)/2),100:end-100,100:end-100,nn));
        im2 = squeeze( v(round(size(v,1)/2),100:end-100,100:end-100,nn + 1));
    elseif strcmp( regdir, 'y' )
        im1 = squeeze( v(100:end-100,round(size(v,2)/2),100:end-100,nn));
        im2 = squeeze( v(100:end-100,round(size(v,2)/2),100:end-100,nn + 1));
    end
    out = ImageCorrelation(im1,im2);
    shift(nn+1) = round( out.shift2);
end

fprintf( '\n Cropping volumes. ' )
z0 = abs( cumsum( shift ) );
z1 = max(z0) - z0;
v2 = zeros( [size(v,1), size(v,2), size(v,3) - max(z0), nn_max], 'uint8' );

for nn = 1:nn_max
    v2(:,:,:,nn) = v(:,:,1+z0(nn):end-z1(nn),nn);
    %disp( size( v(:,:,1+z0(nn):end-z1(nn),nn) ) )
    
end
fprintf( 'Done in %.1f s', toc - t);

%% Save
%nimplay(cat(3,im1(s:end,:),im2(1:end-s+1,:)))
%nimplay( squeeze(v2(xx,:,:,:)));

% Slice: Animated gif
t = toc;
fprintf( '\n Save animated gifs' )

xx = round( size( v2, 1) / 2);

filename = sprintf( '%s/%s_loadSequ_x%04u.gif', scan_path, scan_name, xx );
fprintf( '\n output file: %s', filename)

map = colormap(gray);
for nn = 1:size(v2,4)    
    %[A, map] = gray2ind(rot90(squeeze( v2(xx,:,:,nn))));
    A = gray2ind(rot90(squeeze( v2(xx,:,:,nn))));
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

yy = round( size( v2, 2) / 2);

filename = sprintf( '%s/%s_loadSequ_y%04u.gif', scan_path, scan_name, yy );
fprintf( '\n output file: %s', filename)

for nn = 1:size(v2,4)    
    A = gray2ind(rot90(squeeze( v2(:,yy,:,nn))));
    if nn == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

zz = round( size( v2, 3) / 2);

filename = sprintf( '%s/%s_loadSequ_z%04u.gif', scan_path, scan_name, zz );
fprintf( '\n output file: %s', filename)

for nn = 1:size(v2,4)    
    A = gray2ind((squeeze( v2(:,:,zz,nn))));
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
% save( filename, 'v2' , '-v7.3', '-nocompression')
% fprintf( '\n Done in %.1f s', toc - t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
