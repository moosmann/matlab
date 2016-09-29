function p05_reco(data_dir, bin, poolsize, rotation, subfolder, write_proj, write_reco, verbose)
% P05 reconstruction pipeline.
%
% Written by Julian Moosmann.
% First version: 2016-09-28. Last modifcation: 2016-09-29.
%

%% TODO: ROTCENTER
%% TODO: ROI reconstruction
%% TODO: pixel filtering threshold
%% TODO: make 'read file names into matrix' a function
%% TODO: make read filenames into matrix/struct a function

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if nargin < 1
    %data_dir = pwd;
    data_dir = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160920_000_diana/raw/Mg-10Gd39_1w';
end
if nargin < 2
    bin = 1;
end
if nargin < 3
    poolsize = 30; 
end
if nargin < 4
    rotation = 180;
end
if nargin < 5
    subfolder = 'test';
end
if nargin < 6
    write_proj = 0;
end
if nargin < 7
    write_reco = 1;
end
if nargin < 8
    verbose = 1;
end
%% Hard coded arguments. TODO: make optional
ring_filter = 0;

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PrintVerbose(verbose, '\nStart P05 reconstruction pipeline')
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';

% Extract folder and scan names
while data_dir(end) == filesep
    data_dir = data_dir(1:end-1);
end
[raw_dir, scan] = fileparts(data_dir);
data_dir = [data_dir, filesep];
[beamtime_dir, raw_folder] = fileparts(raw_dir);
if ~strcmp(raw_folder, 'raw')
    fprintf('\n ERROR. Folder name is not raw: %s\n', raw_folder)
    return
end
PrintVerbose(verbose, '\n raw folder: %s', raw_dir)
PrintVerbose(verbose, '\n scan: %s', scan)

% output folder
processed_dir = [beamtime_dir, filesep, 'processed', filesep, scan];
if isempty(subfolder)
    out_dir = processed_dir;
else
    out_dir = [processed_dir, filesep, subfolder];
end
PrintVerbose(verbose, '\n output directory: %s', out_dir)
flatcor_dir = [out_dir, filesep, 'flat_corrected', filesep];
reco_dir = [out_dir, filesep, 'reco', filesep];
CheckAndMakePath(reco_dir) 
CheckAndMakePath(flatcor_dir)

%% File names
% Raw projection file names
data_struct = dir( [data_dir, '*.img'] );
img_names = {data_struct.name};
num_img = numel(img_names);

% Ref file names
data_struct = dir( [data_dir, '*.ref'] );
ref_names = {data_struct.name};
num_ref = numel(ref_names);

% Dark file names
data_struct = dir( [data_dir, '*.dar'] );
dark_names = {data_struct.name};
num_dark = numel(dark_names);

PrintVerbose(verbose, '\n #[dark, ref, img] = [%g, %g, %g]', num_dark, num_ref, num_img)

%% Dark field
t = toc;
PrintVerbose(verbose, '\n Processing dark fields.')
for nn = num_dark:-1:1
    stack(:, :, nn) = FilterPixel( read_dat( sprintf('%s%s', data_dir, dark_names{nn}) ), [0.02 0.02]);
end
dark = FilterPixel( squeeze( median(stack, 3) ), [0.02 0.02]);
[dim1, dim2] = size(dark);
if sum(dark <= 0)
    fprintf('\n CAUTION: dark field contains zeros')
end
PrintVerbose(verbose, ' Elapsed time: %g s', toc-t)
PrintVerbose(verbose, '\n image size: [%g, %g]', dim1, dim2)

%% Flat field
t = toc;
PrintVerbose(verbose, '\n Processing flat fields.')
% Preallocation
stack = zeros(dim1, dim2, num_ref, 'single');
ref_names_mat = NameCellToMat(ref_names);
% Parallel loop
CheckParpool(poolsize);
parfor nn = 1:num_ref    
    stack(:, :, nn) = FilterPixel( read_dat( sprintf('%s%s', data_dir, ref_names_mat(nn, :)) ), [0.01 0.005]);
end
flat = FilterPixel( squeeze( mean(stack, 3) ), [0.005 0.0025]);
% Binning
if bin
    flat = Binning( flat );
    dark = Binning( dark );
end
% Write flat field
if write_proj
    filename = sprintf('%sflat.tif', flatcor_dir);
    write32bitTIFfromSingle(filename, flat);    
end
if sum(flat <= 0)
    fprintf('\n CAUTION: flat field contains zeros')
end
% Subtract dark from flat in place
flat = flat - dark;
if sum(flat <= 0)
    fprintf('\n CAUTION: dark-field-corrected flat field contains zeros')
end
PrintVerbose(verbose, ' Elapsed time: %g s', toc-t)

%% Projections
t = toc;
PrintVerbose(verbose, '\n Processing raw projections.')
stack = zeros( size(flat,1), size(flat,2), num_img, 'single');
img_names_mat = NameCellToMat( img_names );
parfor nn = 1:num_img
    img = FilterPixel( read_dat( sprintf('%s%s', data_dir, img_names_mat(nn, :)) ), [0.02 0.01]);
    if bin
        img = Binning( img );
    end
    img = img - dark;
    img = img ./ flat;
    stack(:, :, nn) = img;
end
fprintf(' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)

%% Ring artifact filter
if ring_filter
    t = toc;
    PrintVerbose(verbose, '\n Ring artifact filtering.')
    sino_mean = mean( stack, 3);
    sino_mean_med = medfilt2( sino_mean, [9, 1], 'symmetric' );
    mask = sino_mean_med ./ sino_mean;
    parfor nn = 1:num_img
        stack(:, :, nn) = stack(:, :, nn) .* mask;
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
end

%% Write corrected projections
if write_proj
    t = toc;
    PrintVerbose(verbose, '\n Saving corrected projections.')
    % Dark
    filename = sprintf('%sdark.tif', flatcor_dir);
    write32bitTIFfromSingle(filename, dark);    
    % Projections
    parfor nn = 1:num_img        
        filename = sprintf('%sproj_%06u.tif', flatcor_dir, nn );
        write32bitTIFfromSingle(filename, stack(:, :, nn) );
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
end

%% Tomgraphic reconstruction and write slices
t = toc;
PrintVerbose(verbose, '\n Tomographic reconstuction.')
theta = rotation * (0:num_img - 1) / (num_img - 1);
output_size = size( stack,1);
rot_axis_pos_offset = 3;
rot_axis_pos = round(size( stack, 1) / 2) + rot_axis_pos_offset;
filt = iradonDesignFilter('Ram-Lak', size(stack, 1), 1 , 0);

for nn = round( size( stack, 2))
    sino = squeeze( stack(:, nn, :) );
    sino_ft = fft( sino, [], 1);
    sino_ft = bsxfun(@times, sino_ft, filt);
    sino = ifft( sino_ft , [], 1, 'symmetric');
       
end

% parfor nn = 1:size( stack, 2)
%     sino = RotAxisSymmetricCropping(squeeze( stack(:, nn, :) ), rot_axis_pos, 1);
%     rec = iradon(sino, theta, 'linear', 'Ram-Lak', 1, output_size);
%     if write_reco
%         filename = sprintf('%sreco_%06u.tif', reco_dir, nn );
%         write32bitTIF(filename, rec)
%     end    
% end


PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PrintVerbose(verbose, '\n Finished. Total time elapsed:%g s = %g min\n', toc-t, (toc-t)/60);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckParpool(poolsize)
%% Parpool
% check if more than 1 worker should be used
if poolsize > 1    
    % check if parpool exists
    if isempty(gcp('nocreate'))
        parpool( poolsize);
    else
        % check number of worker of existing pool
        poolobj = gcp;
        if poolobj.NumWorkers ~= poolsize
            poolobj.delete;
            parpool( poolsize);
        end        
    end
end
end
