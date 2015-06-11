function im_stack=readcellstack(folder,part_of_nameprefix,skip_image,step_size,num_steps,namepostfix,median_filtering)
% Script for reading stack of images (ONLY .tif) from Photron camera during
% phase-stepping technique

% Written by J.Moosmann and V.Altapova, some days ago :) 20-24 of August
% Welcome to modify if you know what to do

% folder - string: path to the tif files, last slash is optional
% part_of_nameprefix - string: optional string. necessary if there are
% tif-files in the directory 'folder' which should not be read and which
% have a different name pattern.  if not necessary choose empty string ('').
% skip_image - integer: how many images to skip in the begining
% step_size - integer: how many images per position of the grating
% namepostfix - string: format of the images. up to now only tif files.
% median_filtering - boolean: simle median filter: takes median 3x3 pixels (decrease the influense of noise)
    
    if nargin<2, part_of_nameprefix = ''; end; 
    if nargin<3, skip_image=1; end;
    if nargin<4, step_size = 1; end;
    if nargin<5, num_steps = 0; end;
    if nargin<6, namepostfix = 'tif'; end;
    if nargin<7, median_filtering = 0; end;
   
% Check if last character of  string 'folder' is a slash (/).
if (folder(end)~='/'), folder = [folder,'/'];end

% Get struct 'files' containing the folder names.
files = dir([folder part_of_nameprefix '*.' namepostfix]);
num_of_files = numel(files);

% Read images into image stack.
% Preallocation of memory improves performance dramatically.
im_stack      = cell(1,num_of_files);
% Define which images to read.
if num_steps==0,
step_list     = skip_image:step_size:num_of_files;
else,
step_list     = skip_image:step_size:skip_image+(num_steps-1)*step_size;
end;
num_of_images = length(step_list);
% Loop over file names files(i).name.
counter = 1;
tic;
for k = step_list,
    % Mandatory conversion to double precision reduces performance.
    if median_filtering==1,
    im_stack{k} = medfilt2(double(imread([folder files(k).name])));    
    else,
    im_stack{k} = double(imread([folder files(k).name]));    
    end;
    counter = counter + 1;
end;
read_time = toc;

% Print paramters.
fprintf(['Read %i of %i ''' namepostfix ''' files in %gs contained in:\n' folder '\n'],num_of_images,num_of_files,read_time);
