% write32bitTIF.m
% writing 32bit float images into a TIFF container
% readability with imagej/fiji is guaranteed, however, 32bit TIFFs seem not
% to be a standard format (imagej/fiji works fine with it)!
% usage:
%           write32bitTIF(filename,A)
%
% Author: Daniel Hänschke
% 06/12/2012, 10:00
% - compatibility with matlab2012a checked

function write32bitTIFfromSingle(filename,A, fmode)
if nargin < 3
    fmode = 'w';
end

% write the image
tiffOut = Tiff(filename, fmode);
tiffOut.setTag('Photometric',1);
tiffOut.setTag('ImageLength',size(A,1));
tiffOut.setTag('ImageWidth',size(A,2));
tiffOut.setTag('BitsPerSample',32);
tiffOut.setTag('SamplesPerPixel',1);
tiffOut.setTag('RowsPerStrip',size(A,1));
tiffOut.setTag('PlanarConfiguration',1);
tiffOut.setTag('SampleFormat',3);
tiffOut.write(A);
tiffOut.close();
end