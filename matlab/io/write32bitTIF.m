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

function write32bitTIF(filename,A)
% check input data
% ----------------
% if ~ismatrix(A)
%     error('Input variable ''A'' is not a matrix.\n');
% end
% if ~isreal(A) 
%     error('Input matrix ''A'' cannot be associated with a real valued image.\n');
% end
% % check output file name
% % ----------------------
% % does the file already exist?
% if exist(filename,'file')
%     error('%s already exists.\n',filename);
% end
% % does the output directory exist?
% parseTmp = strfind(filename,'/');
% if ~isempty(parseTmp)
%     dirCheck = filename(1:parseTmp(end)-1);
%     if ~exist(dirCheck,'dir')
%         error('directory %s doesn''t exist.\n', dirCheck);
%     end    
% end
% write the image
tiffOut = Tiff(filename, 'w');
tiffOut.setTag('Photometric',1);
tiffOut.setTag('ImageLength',size(A,1));
tiffOut.setTag('ImageWidth',size(A,2));
tiffOut.setTag('BitsPerSample',32);
tiffOut.setTag('SamplesPerPixel',1);
tiffOut.setTag('RowsPerStrip',size(A,1));
tiffOut.setTag('PlanarConfiguration',1);
tiffOut.setTag('SampleFormat',3);
tiffOut.write(single(A));
tiffOut.close();
end