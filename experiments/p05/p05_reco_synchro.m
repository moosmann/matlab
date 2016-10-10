% P05 reconstruction pipeline.
%
% Written by Julian Moosmann.
% First version: 2016-09-28. Last modifcation: 2016-10-10

%% TODO: ROT CENTER
%% TODO: ROI reconstruction
%% TODO: pixel filtering threshold
%% TODO: make read filenames into matrix/struct a function
%% TODO: optimize backprojection filter, it's recalculated each time
%% TODO: improve syntax of filter for FBP
%% TODO: add option to save phase maps

%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data_dir = pwd;
data_dir = ...
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_030';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_1000';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_200';    
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_030';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_003_synchro';
proj_stride = 1; % Stride of projection images to read
bin = 1; % bin size: if 2 do 2 x 2 binning, if 1 do nothing
poolsize = 30; % number of workers in parallel pool to be used
gpu_ind = 1; % GPU Device to use: gpuDevice(gpu_ind)
gpu_thresh = 0.7; % Percentage of available GPU memory to be used at maximum
angles = 2 * pi; % full angular range (first to very last image) or array of angles in radians
correlate_proj_and_flat = 1;% correlate flat fields and projection to correct beam shaking
correlate_proj_and_flat_use_best_match = 0; % if 0 takes the mean of all flats which are closest to the projection within the order of 1 pixel
write_proj = 1;
write_reco = 1;
write_to_scratch = 1; % for testing
verbose = 1; % print information to standard output
phase_retrieval_method = ''; % Phase retrieval if phase_retrieval_method not empty
ring_filter = 1; % ring artifact filter
ring_filt_med_width = 11;
vol_size = 0; % number of voxels of the volume to be reconstructed
rotation_axis_offset = []; % if empty use automatic computation
roi1 = [0.25 0.75];
roi2 = [0.25 0.75];
filter_pad_method = 'symmetric';'none';  'none'; % FBP filter padding type
filter_pad_length = 'twice'; % FBP filter padding length
take_log = 0; % logarithm for attenuation contrast
subfolder = 'test_nobin_ringFilter' ; % folder to subfolder flat_corrected, reco, if not empty

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
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
if write_to_scratch
    folder_main = 'scratch_cc';
else
    folder_main = 'processed';
end
processed_dir = [beamtime_dir, filesep, folder_main, filesep, scan];
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
binned_size = floor( raw_size / bin );
PrintVerbose(verbose, '\n raw image size: [%g, %g]', raw_size)
PrintVerbose(verbose, '\n binned image size: [%g, %g]', binned_size)
PrintVerbose(verbose, '\n Processing dark fields.')
dark = zeros( [binned_size, num_dark], 'single');
for nn = 1:num_dark
    dark(:, :, nn) = Binning( FilterPixel( read_dat( sprintf('%s%s', data_dir, dark_names{nn}) ), [0.02 0.02]), bin);    
end
dark = FilterPixel( squeeze( median(dark, 3) ), [0.02 0.02]);
if sum(dark <= 0)
    fprintf('\n CAUTION: dark field contains zeros')
end
PrintVerbose(verbose, ' Elapsed time: %g s', toc-t)

%% ROI
if numel(roi1) == 2
    roi1 = floor( roi1(1) * binned_size(1) ):ceil( roi1(2) * binned_size(1) );
end
if numel(roi2) == 2
    roi2 = floor( roi2(1) * binned_size(2) ):ceil( roi2(2) * binned_size(2) );
end

%% Parallel CPU pool
t = toc;
PrintVerbose( poolsize > 1, '\n Use parallel CPU pool of %u workers (or more). ', poolsize)
OpenParpool(poolsize);
PrintVerbose( poolsize > 1, ' Elapsed time: %g s', toc-t)

%% Flat field
t = toc;
PrintVerbose(verbose, '\n Processing flat fields.')
% Preallocation
flat = zeros( [binned_size, num_ref], 'single');
ref_names_mat = NameCellToMat(ref_names);
% Parallel loop
parfor nn = 1:num_ref
    flat(:, :, nn) = Binning( FilterPixel( read_dat( sprintf('%s%s', data_dir, ref_names_mat(nn, :)) ), [0.01 0.005]), bin) - dark;    
end
flat_mean = FilterPixel( squeeze( mean(flat, 3) ), [0.005 0.0025]);
% Write mean flat field
if write_proj    
    write32bitTIFfromSingle( sprintf('%sflat_mean.tif', flatcor_dir), flat_mean);    
end
if sum(flat_mean <= 0)
    fprintf('\n CAUTION: mean flat field contains zeros')
end
PrintVerbose(verbose, ' Elapsed time: %g s', toc-t)

%% Projections
t = toc;
PrintVerbose(verbose, '\n Read and filter raws.')
% Preallocation
proj = zeros( binned_size(1), binned_size(2), num_proj, 'single');
img_names_mat = NameCellToMat( img_names );
img_names_mat = img_names_mat(proj_ind, :);
parfor nn = 1:num_proj    
    img = Binning( FilterPixel( read_dat( sprintf('%s%s', data_dir, img_names_mat(nn, :)) ), [0.02 0.01]), bin) - dark;
    if ~correlate_proj_and_flat
        img = img ./ flat_mean;
    end
    proj(:, :, nn) = img;
end
PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc - t, ( toc - t ) / 60 )

%% Correlate shifted flat fields
if correlate_proj_and_flat    
    PrintVerbose(verbose, '\n Correct beam shake.')
    flatroi = flat(1:50, roi2, :);
    projroi = proj(1:50, roi2, :);
    xshift = zeros( num_proj, num_ref);
    yshift = xshift;
    for ff = 1:num_ref
        flat_ff = flatroi(:,:,ff);
        parfor pp = 1:num_proj
            out = ImageCorrelation(projroi(:,:,pp), flat_ff, 0,0);
            xshift(pp, ff) = out.Xshift;
            yshift(pp, ff) = out.Yshift; % relevant shift
        end    
    end
    if correlate_proj_and_flat_use_best_match        
        parfor nn = 1:num_proj
            [val, pos] = min( abs( yshift( nn, : ) ) );       
            proj(:, :, nn) = proj(:, :, nn) ./ flat(:, :, pos);
        end
    else
        flat_count = zeros(1, num_proj);
        parfor nn = 1:num_proj
            vec = 1:num_ref;
            flat_ind = [];
            mm = 0;
            while isempty( flat_ind ) && mm < 100; % arbitrary maxim um shift
                flat_ind = vec( abs( round( yshift(nn, :) ) ) == mm );
                mm = mm + 1;       
            end
            flat_count(nn) = numel(flat_ind);    
            flat_selected_mean = FilterPixel( squeeze( mean( flat(:, :, flat_ind), 3) ), [0.005 0.0025]);
            proj(:, :, nn) = proj(:, :, nn) ./ flat_selected_mean;
        end
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc - t, ( toc - t ) / 60 )
    if ~correlate_proj_and_flat_use_best_match
        PrintVerbose(verbose, '\n Number of flats used: [mean, min, max] = [%g, %g, %g]', mean( flat_count ), min( flat_count ), max( flat_count) )
    end
end
PrintVerbose(verbose, '\n sinogram size = [%g, %g, %g]', size( proj ) )

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
        write32bitTIFfromSingle(filename, proj(:, :, nn) );
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
end

%% Ring artifact filter
if ring_filter
    t = toc;
    PrintVerbose(verbose, '\n Ring artifact filtering.')
    sino_mean = mean( proj, 3);
    sino_mean_med = medfilt2( sino_mean, [ring_filt_med_width, 1], 'symmetric' );
    mask = sino_mean_med ./ sino_mean;
    parfor nn = 1:num_proj
        proj(:, :, nn) = proj(:, :, nn) .* mask;
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
end


%% Tomgraphic reconstruction 
filter_type = 'Ram-Lak';
%% TODO: rename
pixel_size = [1, 1];
filter_direction = 1;
t = toc;
if numel( angles ) == 1
    angles = angles * (0:num_img - 1) / (num_img - 1);
end

%% Rotation axis position
% retrieve index at angles 0 and pi
[val1, ind1] = min( abs( angles ));
im1 = proj( :, roi2, ind1);
[val2, ind2] = min( abs( angles - pi ));
im2 = flipud( proj( :, roi2, ind2) );
PrintVerbose(verbose, '\n Rotation axis:\n  Correlate image #%g at angle %g and #%g at angle %g', ind1, val1, ind2, val2 + pi)
% Correlated images
out = ImageCorrelation( im1, im2, 0, 0);
PrintVerbose(verbose, '\n  respective shift: %g, %g', out.Xshift, out.Yshift)
PrintVerbose(verbose, '\n  calulated rotation axis offset: %g, %g', out.Xshift / 2)
if isempty(rotation_axis_offset)
    rotation_axis_offset = round( out.Xshift / 2, 1);
end

%% Phase retrieval
PrintVerbose(verbose, '\n Phase retrieval.')
if ~isempty( phase_retrieval_method )
    pf = PhaseFilter(phase_retrieval_method, binned_size, [30 1 1e-6], 1.5); 
    parfor nn = 1:num_proj
        proj(:, :, nn) = real( ifft2( pf .* fft2( squeeze( proj(:, :, nn) ) ) ) )
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
else
    PrintVerbose(verbose, ' None.')
end

% Determine maximum slab sizes for GPU reconstruction
num_pix1 = size( proj, 1);
num_pix2 = size( proj, 2);
if vol_size == 0
    vol_size = [num_pix1, num_pix1, num_pix2];
end
num_proj =  size( proj, 3);
% GPU memory required for reconstruction
rec_mem = ( num_pix1 * num_proj + vol_size(1) * vol_size(2) ) * 4;
% GPU memory limit
gpu = gpuDevice( gpu_ind );
num_sli = floor( gpu_thresh * gpu.AvailableMemory / rec_mem );
% slabs required
num_slabs = ceil( vol_size(3) / num_sli );
% Readjust number slices to make all of similar size
num_sli = ceil(vol_size(3) / num_slabs);
subvol_size = [vol_size(1:2), num_sli];
PrintVerbose(verbose, '\n rotation axis offset = %g', rotation_axis_offset );
PrintVerbose(verbose, '\n GPU name: %s', gpu.Name );
PrintVerbose(verbose, '\n GPU memory: %g MiB (available), %g MiB (total) (%.2f%%)', gpu.AvailableMemory/10^6, gpu.TotalMemory/10^6, 100*gpu.AvailableMemory/gpu.TotalMemory)
PrintVerbose(verbose, '\n size of total volume to be reconstructed = [%g, %g, %g]', vol_size )
PrintVerbose(verbose, '\n maximum size of subvolume to be reconstructed = [%g, %g, %g]', subvol_size )
PrintVerbose(verbose, '\n estimate of required GPU memory / slice = %g MiB', rec_mem / 1024^2 )
PrintVerbose(verbose, '\n memory of total volume to be reconstructed = %g MiB', prod( vol_size ) * 4 / 1024^2 )
PrintVerbose(verbose, '\n maximum memory of subvolume to be reconstructed = %g MiB', prod( vol_size ) * 4 / 1024^2 )
PrintVerbose(verbose, '\n number of subvolume slabs = %g', num_slabs )
PrintVerbose(verbose, '\n maximum number of slices per slab = %g', num_sli )
PrintVerbose(verbose, '\n Tomographic reconstuction.')



% Loop over slabs
for nn = 1:num_slabs
    
    % Slice indices
    sli0 = 1 + ( nn - 1) * num_sli;
    sli1 = min( [nn * num_sli, vol_size(3)] );
    
    % Filter sinogram
    PrintVerbose(verbose, ' \n Subvolume #%g: filtering sino,', nn)
    sino = FilterSinoForBackproj( permute( proj(:, sli0:sli1, :), [1 3 2] ), filter_direction, filter_type, filter_pad_method, filter_pad_length);
    if take_log
        sino = - reallog( sino );
    end
    
    % Backprojection
    PrintVerbose(verbose, ' backprojection,')
    astra_clear
    subvol_size(3) = size( sino, 3);
    vol = astra_parallel3D( sino, angles, rotation_axis_offset, subvol_size, pixel_size);
    
    % Save subvolume
    if write_reco        
        PrintVerbose(verbose, ' saving.')
        parfor ii = 1:size( vol, 3)
            filename = sprintf( '%sreco_%06u.tif', reco_dir, sli0 + ii - 1);
            write32bitTIF( filename, vol( :, :, ii) )
        end
    end
end
PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60 )

PrintVerbose(verbose, '\n Finished. Total time elapsed: %g s = %g min\n', toc, toc / 60 );
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
