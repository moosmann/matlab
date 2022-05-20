% Skript converting images to dicom format

p = '/asap3/petra3/gpfs/p05/2021/data/11008741/processed/syn0114_72R_Mg10Gd_12w_000/reco/';
d = dir( [p 'float_rawBin4/*tif'] );
fprintf( '\n folder: %s', p )
fprintf( '\n %u tif found', numel(d) )
fprintf( '\n' )

volume_min = -0.0268398;
volume_max = 0.0426215;

for n = 200:760
    
    imname = [d(n).folder filesep d(n).name];
    im = imread(imname);
%     domain(im)
    im = (im - volume_min) ./ (volume_max - volume_min);
    domain(im)
    im = uint16((2^16 - 1) * im );
%     domain(im)
    tmp = d(n).name;    
    dicomname = [p 'dicom/' tmp(1:end-3) 'dcm'];
%     disp(imname);
%     disp(dicomname);
    status = dicomwrite( im, dicomname);
end