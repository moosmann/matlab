function im = read_image(image_path,filetype)
% Reads 'tif' and 'edf' images and all the files MATLAB's imread function
% can read without additional arguments than file format
    
%image_path - string: absolute filename
%filetype - string: namepostfix of image file
    

if strcmpi(filetype,'edf')
    im = pmedfread(image_path);
else
    im = double(imread(image_path,filetype));
end;
