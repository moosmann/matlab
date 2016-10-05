% P05 reconstruction pipeline.
%
% Written by Julian Moosmann.
% First version: 2016-09-28. Last modifcation: 2016-10-05.

%% TODO: ROT CENTER
%% TODO: ROI reconstruction
%% TODO: pixel filtering threshold
%% TODO: make 'read file names into matrix' a function
%% TODO: make read filenames into matrix/struct a function

tic
%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data_dir = pwd;
data_dir = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160920_000_diana/raw/Mg-10Gd39_1w';
proj_stride = 1; % Stride of projection images to read
bin = 0; % if 1, do 2 x 2 binning, else nothing
poolsize = 30;
gpu_ind = 1; % GPU Device to use: gpuDevice(gpu_ind)
gpu_thresh = 0.5; % Percentage of available GPU memory to be used at maximum
angles = pi;
subfolder = 'test';
write_proj = 0;
write_reco = 0;
verbose = 1;
ring_filter = 0;
vol_size = 0; % number of voxels of the volume to be reconstructed
rotation_axis_offset = 3;

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrintVerbose(verbose, '\nStart reconstruction pipeline for P05 ')
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';

%% Input folder
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

%% Output folder
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

% Projection to read
proj_stride = max( 1, proj_stride);
proj_ind = 1:proj_stride:num_img;
num_proj = numel( proj_ind );
PrintVerbose(verbose, '\n Projections: stride = %g, #proj = %g', proj_stride, num_proj)

%% Dark field
t = toc;
raw_size = size( read_dat( sprintf('%s%s', data_dir, dark_names{1}) ) );
PrintVerbose(verbose, '\n raw image size: [%g, %g]', raw_size)
PrintVerbose(verbose, '\n Processing dark fields.')
stack = zeros( [raw_size, num_dark], 'single');
for nn = 1:num_dark
    stack(:, :, nn) = FilterPixel( read_dat( sprintf('%s%s', data_dir, dark_names{nn}) ), [0.02 0.02]);
end
dark = FilterPixel( squeeze( median(stack, 3) ), [0.02 0.02]);
if sum(dark <= 0)
    fprintf('\n CAUTION: dark field contains zeros')
end
PrintVerbose(verbose, ' Elapsed time: %g s', toc-t)

%% Parallel CPU pool
t = toc;
PrintVerbose( poolsize > 1, '\n Use parallel CPU pool of %u workers (or more). ', poolsize)
OpenParpool(poolsize);
PrintVerbose( poolsize > 1, ' Elapsed time: %g s', toc-t)

%% Flat field
t = toc;
PrintVerbose(verbose, '\n Processing flat fields.')
% Preallocation
stack = zeros( [raw_size, num_ref], 'single');
ref_names_mat = NameCellToMat(ref_names);
% Parallel loop
parfor nn = 1:num_ref    
    stack(:, :, nn) = FilterPixel( read_dat( sprintf('%s%s', data_dir, ref_names_mat(nn, :)) ), [0.01 0.005]);
end
flat = FilterPixel( squeeze( mean(stack, 3) ), [0.005 0.0025]);
% Binning
if bin
    flat = Binning( flat );
    dark = Binning( dark );
    raw_size_binned = size( flat );    
else
    raw_size_binned = raw_size;
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
PrintVerbose(verbose, '\n binned image size: [%g, %g]', raw_size_binned)

%% Projections
t = toc;
PrintVerbose(verbose, '\n Processing raw projections.')
% Preallocation
stack = zeros( size(flat,1), size(flat,2), num_proj, 'single');
img_names_mat = NameCellToMat( img_names );
parfor nn = 1:num_proj    
    img = FilterPixel( read_dat( sprintf('%s%s', data_dir, img_names_mat(proj_ind(nn), :)) ), [0.02 0.01]);
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
if numel( angles ) == 1
    angles = angles * (0:num_img - 1) / (num_img - 1);
end
% Determine maximum slab sizes for GPU reconstruction
num_pix1 = size( stack, 1);
num_pix2 = size( stack, 2);
if vol_size == 0
    vol_size = [num_pix1, num_pix1, num_pix2];
end
num_proj =  size( stack, 3);
% GPU memory required for reconstruction
rec_mem = ( num_pix1 * num_proj + vol_size(1) * vol_size(2) ) * 4;
% GPU memory limit
gpu = gpuDevice( gpu_ind );
num_sli = floor( gpu_thresh * gpu.AvailableMemory / rec_mem );
% slabs required
num_slabs = ceil( vol_size(3) / num_sli );
% Readjust slices to get the last slab similar to the other
num_sli = ceil(vol_size(3) / num_slabs);
PrintVerbose(verbose, '\n GPU name: %s', gpu.Name );
PrintVerbose(verbose, '\n GPU memory: %g mb of %g mb available (%.2f%%)', gpu.AvailableMemory/10^6, gpu.TotalMemory/10^6, 100*gpu.AvailableMemory/gpu.TotalMemory)
PrintVerbose(verbose, '\n Tomographic reconstuction.')
for nn = 1:num_slabs
    
    % Slice indices
    sli0 = 1 + ( nn - 1) * num_sli;
    sli1 = min( [nn * num_sli, vol_size(3)] );
    
    % Filter sinogram
    sino = FilterSinoForBackproj( permute( stack(:, sli0:sli1, :), [1 3 2] ), 1, 'Ram-Lak');
    
    % Backprojection
    vol = astra_parallel3D( permute(sino, [1 2 3]), angles, rotation_axis_offset);
    
%     if write_reco
%         filename = sprintf('%sreco_%06u.tif', reco_dir, nn );
%         write32bitTIF(filename, vol)
%     end
end
PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60 )

PrintVerbose(verbose, '\n Finished. Total time elapsed: %g s = %g min\n', toc, toc / 60 );
