% P05 reconstruction pipeline.
%
% Written by Julian Moosmann.
% First version: 2016-09-28. Last modifcation: 2016-10-07

%% TODO: automatic determination of rot center
%% TODO: ROI volume reconstruction
%% TODO: make pixel filtering thresholds parameters
%% TODO: make read filenames into matrix/struct a function
%% TODO: optimize backprojection filter, it's recalculated each time
%% TODO: improve syntax of filter for FBP
%% TODO: correct for beam current
%% TODO: write log file

%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data_dir = pwd;
data_dir = ...
    '/asap3/petra3/gpfs/p05/2016/data/11001464/raw/pnl_16_petrosia_c';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_030';
%'/asap3/petra3/gpfs/p05/2016/commissioning/c20160920_000_diana/raw/Mg-10Gd39_1w';
proj_stride = 1; % Stride of projection images to read
bin = 2; % bin size: if 2 do 2 x 2 binning, if 1 do nothing
poolsize = 28; % number of workers in parallel pool to be used
gpu_ind = 1; % GPU Device to use: gpuDevice(gpu_ind)
gpu_thresh = 0.7; % Percentage of available GPU memory to be used at maximum
angles = 2* pi;[]; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically
subfolder = 'test_ring_filter'; % folder to subfolder flat_corrected, reco, if not empty
correct_beam_shake = 1;% correlate flat fields and projection to correct beam shaking
correct_beam_shake_max_shift = 0.5; % if 0 the best match (i.e. closest to zero) is used, if > 0 all flats which are shifted less than correct_beam_shake_max_shift are used
roi1 = [0.25 0.75]; % for correlation
roi2 = [0.25 0.75]; % for correlation
write_proj = 1;
write_reco = 1;
write_to_scratch = 0; % for testing
verbose = 1; % print information to standard output
ring_filter = 1; % ring artifact filter
ring_filt_med_width = 11;
do_tomo = 1; % reconstruct volume
vol_shape = []; % shape of the volume to be reconstructed, either in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
rotation_axis_offset = 614;[]; % if empty use automatic computation
filter_pad_method = 'none';'symmetric'; 'replicate';  % FBP filter padding type
filter_pad_length = 'nextpow2'; % FBP filter padding length
link_data = 1; % ASTRA data objects become references to Matlab arrays.
take_log = 0; % logarithm for attenuation contrast
butterworth_cutof_frequ = .5;

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
PrintVerbose(verbose, '\nStart reconstruction pipeline for P05 ')
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';
astra_clear % in the case reco was aborted and data objects in ASTRA memory not deleted

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
PrintVerbose(verbose, '\n scan: %s', scan)
PrintVerbose(verbose, '\n raw folder: %s', raw_dir)

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
PrintVerbose(verbose, '\n projections used: %g, stride: %g', num_proj, proj_stride)

%% Dark field
t = toc;
raw_size = size( read_dat( sprintf('%s%s', data_dir, dark_names{1}) ) );
binned_size = floor( raw_size / bin );
PrintVerbose(verbose, '\n image size: raw = [%g, %g], binned = [%g, %g]', raw_size, binned_size)
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
%roi1 = floor( 0.25 * binned_size(1) ):ceil( 0.75 * binned_size(1) );
%roi2 = floor( 0.25 * binned_size(2) ):ceil( 0.75 * binned_size(2) );

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
% if write_proj
%     CheckAndMakePath( flatcor_dir )
%     write32bitTIFfromSingle( sprintf('%sflat_mean.tif', flatcor_dir), flat_mean);    
% end
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
    if ~correct_beam_shake
        img = img ./ flat_mean;
    end
    proj(:, :, nn) = img;
end
PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc - t, ( toc - t ) / 60 )

%% Correlate shifted flat fields
if correct_beam_shake    
    PrintVerbose(verbose, '\n Correct beam shake.')
    % ROI to correlate
    flatroi = flat(1:50, roi2, :);
    projroi = proj(1:50, roi2, :);
    xshift = zeros( num_proj, num_ref);
    yshift = xshift;
    % Compute shift for each pair projection/flat field
    for ff = 1:num_ref
        flat_ff = flatroi(:,:,ff);
        parfor pp = 1:num_proj
            out = ImageCorrelation(projroi(:,:,pp), flat_ff, 0,0);
            xshift(pp, ff) = out.Xshift;
            yshift(pp, ff) = out.Yshift; % relevant shift
        end    
    end
    % use best match
    if correct_beam_shake_max_shift == 0        
        parfor nn = 1:num_proj
            [val, pos] = min( abs( yshift(nn, : ) ) );       
            proj(:, :, nn) = proj(:, :, nn) ./ flat(:, :, pos);
        end
    % use all flats which are shifted less than correct_beam_shake_max_shift pixels
    elseif correct_beam_shake_max_shift > 0
        flat_count = zeros(1, num_proj);
        parfor nn = 1:num_proj
            vec = 1:num_ref;
            flat_ind = vec( abs( yshift(nn, :) ) < correct_beam_shake_max_shift );
            if isempty( flat_ind )
                [~, flat_ind] = min( abs( yshift(nn, :) ) );
            end
            flat_count(nn) = numel(flat_ind);    
            flat_selected_mean = FilterPixel( squeeze( mean( flat(:, :, flat_ind), 3) ), [0.005 0.0025]);
            proj(:, :, nn) = proj(:, :, nn) ./ flat_selected_mean;
        end
    else
        error('Value of maximum shift (%g) is not >= 0', correct_beam_shake_max_shift)
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc - t, ( toc - t ) / 60 )
    if correct_beam_shake_max_shift > 0
        PrintVerbose(verbose, '\n number of flats used per projection: [mean, min, max] = [%g, %g, %g]', mean( flat_count ), min( flat_count ), max( flat_count) )
    end
end
PrintVerbose(verbose, '\n sinogram size = [%g, %g, %g]', size( proj ) )

%% Ring artifact filter
if ring_filter
    t = toc;
    PrintVerbose(verbose, '\n Ring artifact filtering.')
    sino_mean = mean( proj, 3);
    sino_mean_med = medfilt2( sino_mean, [ring_filt_med_width, 1], 'symmetric' );
    mask = sino_mean_med ./ sino_mean;      
    proj = bsxfun( @times, proj, mask);            
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
end

%% Write corrected projections
if write_proj
    t = toc;   
    PrintVerbose(verbose, '\n Saving corrected projections.')
    CheckAndMakePath( flatcor_dir )
    % Dark
%     filename = sprintf('%sdark.tif', flatcor_dir);
%     write32bitTIFfromSingle(filename, dark);    
    % Projections
    parfor nn = 1:num_proj     
        filename = sprintf('%sproj_%06u.tif', flatcor_dir, nn );
        write32bitTIFfromSingle(filename, proj(:, :, nn)' );
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
end

if do_tomo    
   
    %% Tomgraphic reconstruction 
    filter_type = 'Ram-Lak';
    %% TODO: rename
    pixel_size = [1, 1];
    filter_direction = 1;
    t = toc;

    %% Rotation axis position
    PrintVerbose(verbose, '\n Rotation axis:')
    % Automatic determination of full rotation angle if angles is empty
    if isempty( angles )    
        im1 = proj( roi1, roi2, 1);
        im2 = proj( roi1, roi2, num_proj);
        [~, cm1] = ImageCorrelation( im1, im2); % large if full angle is  2 * pi
        [~, cm2] = ImageCorrelation( im1, flipud(im2)); % large if full angle is pi
        if max( cm1(:) ) > max( cm2(:) )
            angles = 2 * pi;
        elseif max( cm1(:) ) < max( cm2(:) )
            angles = pi;
        else
            error('Automatic determination of full angle of rotation via correlation of first and last (flipped) projection failled.')
        end
    end
    if numel( angles ) == 1
        PrintVerbose(verbose, '\n  angular range: %g', angles)
        angles = angles * (0:num_proj - 1) / (num_proj - 1);
    else
        PrintVerbose(verbose, '\n  angular range:', max( angles ) - min( angles ))
    end
    if numel( angles ) ~= num_proj
        error('Number of elements in array of angles (%g) unequal number of projections read (%g)', numel( angles ), num_proj)
    end
    % retrieve index at angles 0 and pi
    [val1, ind1] = min( abs( angles ));
    [val2, ind2] = min( abs( angles - pi ));
    % Correlate images
    im1 = proj( :, roi2, ind1);
    im2 = flipud( proj( :, roi2, ind2) );
    PrintVerbose(verbose, '\n  image to correlate: #%g at angle %g pi and #%g at angle %g pi', ind1, val1 / pi, ind2, (val2 + pi) / pi )
    out = ImageCorrelation( im1, im2, 0, 0);
    PrintVerbose(verbose, '\n  respective shift: %g, %g', out.Xshift, out.Yshift)
    PrintVerbose(verbose, '\n  calulated rotation axis offset: %g, %g', out.Xshift / 2)
    if isempty(rotation_axis_offset)
        rotation_axis_offset = round( out.Xshift / 2, 1);
    end

    %% GPU slab sizes
    num_pix1 = size( proj, 1);
    num_pix2 = size( proj, 2);
    if isempty( vol_shape )
        % default volume given by the detector width and height
        vol_shape = [num_pix1, num_pix1, num_pix2];
    else            
        if vol_shape(1) <  1
            vol_shape(1) = vol_shape(1) * num_pix1;
        end
        if vol_shape(2) <  1
            vol_shape(2) = vol_shape(2) * num_pix1;
        end
        if vol_shape(3) <  1
            vol_shape(3) = vol_shape(3) * num_pix2;
        end
    end
    num_proj =  size( proj, 3);
    % GPU memory required for reconstruction
    rec_mem = ( num_pix1 * num_proj + vol_shape(1) * vol_shape(2) ) * 4;
    % GPU memory limit
    gpu = gpuDevice( gpu_ind );
    num_sli = floor( gpu_thresh * gpu.AvailableMemory / rec_mem );
    % slabs required
    num_slabs = ceil( vol_shape(3) / num_sli );
    % Readjust number slices to make all of similar size
    num_sli = ceil(vol_shape(3) / num_slabs);
    subvol_size = [vol_shape(1:2), num_sli];
    PrintVerbose(verbose, '\n rotation axis offset = %g', rotation_axis_offset );
    PrintVerbose(verbose, '\n GPU name: %s', gpu.Name );
    PrintVerbose(verbose, '\n GPU memory available: %g MiB (%.2f%%) of %g MiB', gpu.AvailableMemory/10^6, 100*gpu.AvailableMemory/gpu.TotalMemory, gpu.TotalMemory/10^6)
    PrintVerbose(verbose, '\n shape of reconstruction volume = [%g, %g, %g]', vol_shape )
    PrintVerbose(verbose, '\n shape of reconstruction subvolume = [%g, %g, %g]', subvol_size )
    PrintVerbose(verbose, '\n estimated GPU memory / slice = %g MiB', rec_mem / 1024^2 )
    PrintVerbose(verbose, '\n memory of total reconstruction volume = %g MiB', prod( vol_shape ) * 4 / 1024^2 )
    PrintVerbose(verbose, '\n maximum memory of subvolume = %g MiB', prod( subvol_size ) * 4 / 1024^2 )
    PrintVerbose(verbose, '\n number of subvolume slabs = %g', num_slabs )
    PrintVerbose(verbose, '\n maximum number of slices per slab = %g', num_sli )
    PrintVerbose(verbose, '\n Tomographic reconstruction.')
    % Loop over slabs
    for nn = 1:num_slabs

        % Slice indices
        sli0 = 1 + ( nn - 1) * num_sli;
        sli1 = min( [nn * num_sli, vol_shape(3)] );

        % Filter sinogram
        PrintVerbose(verbose, ' \n  subvolume %g of %g: Filter sino: ', nn, num_slabs)
        sino = FilterSinoForBackproj( permute( proj(:, sli0:sli1, :), [1 3 2] ), filter_direction, filter_type, filter_pad_method, filter_pad_length, butterworth_cutof_frequ);
        if take_log
            sino = - reallog( sino );
        end
        PrintVerbose(verbose, 'done.')

        % Backprojection
        PrintVerbose(verbose, ' Backprojection: ')    
        subvol_size(3) = size( sino, 3);
        vol = astra_parallel3D( sino, angles, rotation_axis_offset, subvol_size, pixel_size, link_data);
        PrintVerbose(verbose, 'done.')

        % Save subvolume
        if write_reco
            PrintVerbose(verbose, ' Saving:')
            CheckAndMakePath(reco_dir) 
            parfor ii = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', reco_dir, sli0 + ii - 1);
                write32bitTIFfromSingle( filename, vol( :, :, ii) )
            end
            PrintVerbose(verbose, ' done.')
        end
    end
    PrintVerbose(verbose, '\n  Elapsed time: %g s = %g min', toc-t, (toc-t)/60 )

end

PrintVerbose(verbose, '\n Finished. Total time elapsed: %g s = %g min\n', toc, toc / 60 );
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
