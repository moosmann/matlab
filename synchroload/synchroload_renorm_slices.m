input_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/syn166x_3dstacks';
out_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/syn166x_3dstacks_sliceNormalized';

fs = dir( [input_path filesep '*.tif'] );

nn = 1;
filename = [fs(nn).folder filesep fs(nn).name];
tinf = imfinfo( filename );
width = tinf.Width;
height = tinf.Height;
frames = numel( tinf );

CheckAndMakePath( out_path )

vol = zeros( [height width frames], 'double' );
for nn = 1:numel( fs )
    filename = [fs(nn).folder filesep fs(nn).name];
    %% READ
    for mm = 1:frames
        im = imread( filename, mm, 'Info', tinf );
        vol(:,:,mm) = double( im );
    end
    
    %% Normalize
    for ll = height:-1:1
        im = squeeze( vol(ll,:,:) );
        mask = im > 150 & im < 170;
        m(ll) = mean2( im(mask) );
        s(ll) = std( im(mask) );
    end
    m = mean( m ) ./ m;
    s = mean( s ) ./ s;
    
    
    vol = bsxfun( @times, vol, m' );
    
    figure(1), imsc( im ), axis tight
    figure(2), plot( m, '.' ), axis tight
    figure(3), plot( s, '.' ), axis tight
    drawnow
    pause(1)
    
    %% Save
    filename_out = [out_path filesep fs(nn).name];
    im = uint8( vol(:,:,1) );
    imwrite( im, filename_out);
    for mm = 2:frames
        im = uint8( vol(:,:,mm) );
        imwrite( im, filename_out, 'WriteMode','append');        
    end
    
end









whos vol