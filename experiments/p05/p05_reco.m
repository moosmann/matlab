% P05 reconstruction pipeline.
%
% Written by Julian Moosmann.
% First version: 2016-09-28. Last modifcation: 2016-11-08

clear all

%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scan_dir = pwd;
scan_dir = ...      
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_28_00';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1000';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1400';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/mg5gd_21_3w';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_030';
    '/asap3/petra3/gpfs/p05/2016/data/11001464/raw/pnl_16_petrosia_c';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160920_000_diana/raw/Mg-10Gd39_1w';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20161024_000_xeno/raw/xeno_01_b';
read_proj = 0; % Read flatfield-corrected images from disc
read_proj_folder = []; % subfolder of 'flat_corrected' containing projections
proj_stride = 1; % Stride of projection images to read
bin = 4; % bin size: if 2 do 2 x 2 binning, if 1 do nothing
poolsize = 28; % number of workers in parallel pool to be used
gpu_ind = 1; % GPU Device to use: gpuDevice(gpu_ind)
gpu_thresh = 0.8; % Percentage of maximally used to available GPU memory
full_angular_range = [pi]; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically
num_proj_hard = 550; % hard coded number of projection.
correct_beam_shake = 0;%  correlate flat fields and projection to correct beam shaking
correct_beam_shake_max_shift = 0; % if 0: use the best match (i.e. the one which is closest to zero), if > 0 all flats which are shifted less than correct_beam_shake_max_shift are used
write_proj = 0;
write_reco = 1;
write_to_scratch = 1; % write to 'scratch_cc' instead of 'processed'
ring_filter = 1; % ring artifact filter
ring_filt_med_width = 11;
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';
energy = 25; % in keV
sample_detector_distance = 1.0; % in m
phys_pixel_size = bin * 0.68e-6; %1.3056e-5;% physical size of detector pixels, in m
reg_par = 3;
bin_filt = 0.1;
do_tomo = 1; % reconstruct volume
vol_shape = []; % shape of the volume to be reconstructed, either in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
vol_size = []; % if empty, unit voxel size is assumed
rotation_axis_offset = [5.5]; % if empty use automatic computation
rot_global = pi; % global rotation of reconstructed volume
rot_axis_roi1 = [0.25 0.75]; % for correlation
rot_axis_roi2 = [0.25 0.75]; % for correlation
filter_type = 'Ram-Lak';
pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted 
link_data = 1; % ASTRA data objects become references to Matlab arrays.
take_neg_log = 0; % logarithm for attenuation contrast
parfolder = ''; % parent folder to 'reco' and 'flat_corrected'
parfolder_flatcor = ''; % parent folder to 'flat_corrected'
parfolder_reco = ''; % parent folder to 'reco'
verbose = 1; % print information to standard output

%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: automatic determination of rot center (entropy type)
% TODO: make pixel filtering thresholds parameters
% TODO: vertical ROI reco
% TODO: convert read filenames into matrix/struct a function
% TODO: padding options for FBP filter
% TODO: normalize proj with beam current
% TODO: write log file
% TODO: check offsets in projection correlation for rotation axis determination
% TODO: output file format option
% TODO: excentric rotation axis 
% TODO: interactive determination of rotation axis position
% TODO: save and read sinograms
% TODO: read number of projection from log file
% TODO: add photometric tag when reading tif files without one
% TODO: read imager ROI and check if it's faster at all
% TODO: optional save format. Currently 32 bit tiff.
% TODO: ref stride

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
warning( 'off','all')
cdscandir = cd(scan_dir);
PrintVerbose(verbose, '\nStart P05 reconstruction pipeline of scan: ')
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';
astra_clear % if reco was aborted, data objects are not deleted from ASTRA memory

%% Input folder
while scan_dir(end) == filesep    
    scan_dir(end) = [];
end
[raw_dir, scan] = fileparts(scan_dir);
scan_dir = [scan_dir, filesep];
[beamtime_dir, raw_folder] = fileparts(raw_dir);
if ~strcmp(raw_folder, 'raw')
    fprintf('\n ERROR. Folder name is not raw: %s\n', raw_folder)
    return
end
PrintVerbose(verbose, '%s', scan)
PrintVerbose(verbose, '\n scan dir:%s', scan_dir)

%% Output folder
out_folder = 'processed'; 
if write_to_scratch
    out_folder = 'scratch_cc';
end
processed_dir = [beamtime_dir, filesep, out_folder, filesep, scan];
if isempty(parfolder)
    out_dir = processed_dir;
else
    out_dir = [processed_dir, filesep, parfolder];
end
% flat_corrected dir
if isempty( parfolder_flatcor )
    flatcor_dir = [out_dir, filesep, 'flat_corrected', filesep];
else
    flatcor_dir = [out_dir, filesep, 'flat_corrected', filesep, parfolder_flatcor, filesep];
end
PrintVerbose(verbose & write_proj, '\n flatcor: %s', flatcor_dir)
% reco dir
if isempty( parfolder_reco )
    reco_dir = [out_dir, filesep, 'reco', filesep];
else
    reco_dir = [out_dir, filesep, 'reco', filesep, parfolder_reco, filesep];
end
PrintVerbose(verbose, '\n reco   : %s', reco_dir)

%% File names

% Raw projection file names
data_struct = dir( [scan_dir, '*.img'] );
if isempty( data_struct )
    data_struct = dir( [scan_dir, 'proj_*.tif'] );
end
img_names = {data_struct.name};
img_nums = CellString2Vec( img_names );
num_img = numel(img_names);

% Ref file names
data_struct = dir( [scan_dir, '*.ref'] );
if isempty( data_struct )
    data_struct = dir( [scan_dir, 'ref_*.tif'] );
end
ref_names = {data_struct.name};
ref_nums = CellString2Vec( ref_names );
num_ref = numel(ref_names);

% Dark file names
data_struct = dir( [scan_dir, '*.dar'] );
if isempty( data_struct )
    data_struct = dir( [scan_dir, 'dark_*.tif'] );
end
dark_names = {data_struct.name};
dark_nums = CellString2Vec( dark_names );
num_dark = numel(dark_names);

PrintVerbose(verbose, '\n number of [dark, ref, img] = [%g, %g, %g]', num_dark, num_ref, num_img)

% Projection to read
proj_stride = max( 1, proj_stride);
proj_ind = 1:proj_stride:num_img;
num_proj = numel( proj_ind );
PrintVerbose(verbose, '\n projections used: %g, indices: first:stride:last =  %g:%g:%g', num_proj, proj_ind(1), proj_stride, proj_ind(end))

%% Image shape and ROI
filename = sprintf('%s%s', scan_dir, dark_names{1});
im = read_image( filename );
raw_size = size( im );
binned_size = floor( raw_size / bin );
PrintVerbose(verbose, '\n image shape: raw = [%g, %g], binned = [%g, %g]', raw_size, binned_size)
if numel(rot_axis_roi1) == 2
    rot_axis_roi1 = floor( rot_axis_roi1(1) * binned_size(1) ):ceil( rot_axis_roi1(2) * binned_size(1) );
end
if numel(rot_axis_roi2) == 2
    rot_axis_roi2 = floor( rot_axis_roi2(1) * binned_size(2) ):ceil( rot_axis_roi2(2) * binned_size(2) );
end

%% Start parallel CPU pool
t = toc;
PrintVerbose( poolsize > 1, '\n Use parallel CPU pool of %u workers (or more). ', poolsize)
OpenParpool(poolsize);
PrintVerbose( poolsize > 1, ' Elapsed time: %g s', toc-t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read flat corrected projection
if read_proj    
data_struct = dir( [flatcor_dir, 'proj*.*'] );
if isempty( data_struct )
    
    fprintf('\n No flat corrected projections found! Switch to standard pre-processing.')
    read_proj = 0;
else
    proj_names = {data_struct.name};
    num_proj_read = numel(proj_names);
    if num_proj > num_proj_read
        fprintf('\n Less projections available (%g) than demanded (%g)! Switch to standard pre-processing.', num_proj_read, num_proj )
        read_proj = 0;    
    end
end    
    % File names
    t = toc;
    PrintVerbose(verbose, '\n Read flat corrected projections.')    
    if num_proj_read ~= num_proj
        fprintf('\n CAUTION: Number of flat corrected projections read (%g) differs from number of projections to be processed (%g)!\n', num_proj_read, num_proj)
    end

    % Preallocation
    proj = zeros( binned_size(1), binned_size(2), num_proj_read, 'single');
    proj_names_mat = NameCellToMat( proj_names );

    % Read projections    
    parfor nn = 1:num_proj_read
        filename = sprintf('%s%s', flatcor_dir, proj_names_mat(nn, :));
        %proj(:, :, nn) = imread( filename, 'tif' )';
        proj(:, :, nn) = read_image( filename )';
    end    
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc - t, ( toc - t ) / 60 )

%% Read raw data
elseif ~read_proj
    %% Dark field
    t = toc;
    PrintVerbose(verbose, '\nProcessing dark fields.')
    dark = zeros( [binned_size, num_dark], 'single');
    parfor nn = 1:num_dark
        filename = sprintf('%s%s', scan_dir, dark_names{nn});
        dark(:, :, nn) = Binning( FilterPixel( read_image( filename ), [0.02 0.02]), bin);    
    end
    dark = FilterPixel( squeeze( median(dark, 3) ), [0.02 0.02]);
    if sum(dark <= 0)
        fprintf('\n CAUTION: dark field contains zeros')
    end
    PrintVerbose(verbose, ' Elapsed time: %g s', toc-t)

    %% Flat field
    t = toc;
    PrintVerbose(verbose, '\nProcessing flat fields.')
    % Preallocation
    flat = zeros( [binned_size, num_ref], 'single');
    ref_names_mat = NameCellToMat(ref_names);
    % Parallel loop
    parfor nn = 1:num_ref
        filename = sprintf('%s%s', scan_dir, ref_names_mat(nn, :));
        flat(:, :, nn) = Binning( FilterPixel( read_image( filename ), [0.01 0.005]), bin) - dark;    
    end
    flat_mean = FilterPixel( squeeze( mean(flat, 3) ), [0.005 0.0025]);    
    if sum(flat_mean <= 0)
        fprintf('\n CAUTION: mean flat field contains zeros')
    end
    PrintVerbose(verbose, ' Elapsed time: %g s', toc-t)

    %% Projections
    t = toc;
    PrintVerbose(verbose, '\nRead and filter raws.')
    % Preallocation
    proj = zeros( binned_size(1), binned_size(2), num_proj, 'single');
    img_names_mat = NameCellToMat( img_names );
    img_names_mat = img_names_mat(proj_ind, :);
    parfor nn = 1:num_proj
        filename = sprintf('%s%s', scan_dir, img_names_mat(nn, :));
        img = Binning( FilterPixel( read_image( filename ), [0.02 0.01]), bin) - dark;
        if ~correct_beam_shake
            img = img ./ flat_mean;
        end
        proj(:, :, nn) = img;
    end
    PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc - t, ( toc - t ) / 60 )

    %% Correlate shifted flat fields
    if correct_beam_shake    
        PrintVerbose(verbose, '\nCorrect beam shake.')
        % ROI to correlate
        flatroi = flat(1:50, rot_axis_roi2, :);
        projroi = proj(1:50, rot_axis_roi2, :);
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
        % use all flats which are shifted less pixels than correct_beam_shake_max_shift
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
        PrintVerbose(verbose, '\nFilter ring artifacts.')
        sino_mean = mean( proj, 3);
        sino_mean_med = medfilt2( sino_mean, [ring_filt_med_width, 1], 'symmetric' );
        mask = sino_mean_med ./ sino_mean;      
        proj = bsxfun( @times, proj, mask);            
        PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
    end

    %% Write corrected projections
    if write_proj
        t = toc;   
        PrintVerbose(verbose, '\nSave flat-corrected projections.')
        CheckAndMakePath( flatcor_dir )
        % Projections
        parfor nn = 1:num_proj     
            filename = sprintf('%sproj_%06u.tif', flatcor_dir, nn );
            write32bitTIFfromSingle(filename, proj(:, :, nn)' );
        end
        PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  %% Rotation axis position and Tomgraphic parameters %%%%%%%%%%%%%%%%%%%
if do_tomo    
    t = toc;
   
    PrintVerbose(verbose, '\n Rotation axis:')
    % Automatic determination of full rotation angle if full_angular_range is empty
    if isempty( full_angular_range )    
        im1 = proj( rot_axis_roi1, rot_axis_roi2, 1);
        im2 = proj( rot_axis_roi1, rot_axis_roi2, num_proj);
        [~, cm1] = ImageCorrelation( im1, im2); % large if full angle is  2 * pi
        [~, cm2] = ImageCorrelation( im1, flipud(im2)); % large if full angle is pi
        if max( cm1(:) ) > max( cm2(:) )
            full_angular_range = 2 * pi;
        elseif max( cm1(:) ) < max( cm2(:) )
            full_angular_range = pi;
        else
            error('Automatic determination of full angle of rotation via correlation of first and last (flipped) projection failed.')
        end
    end
    %if numel( angles ) == 1
    PrintVerbose(verbose, '\n  angular range: %g', full_angular_range)
    %angles = full_angular_range * (0:num_proj - 1) / (num_proj - 1);
    angles = full_angular_range * img_nums / num_proj_hard;
%     else
%         PrintVerbose(verbose, '\n  angular range:', max( angles ) - min( angles ))
%     end
    if numel( angles ) ~= num_proj
        error('Number of elements in array of angles (%g) unequal number of projections read (%g)', numel( angles ), num_proj)
    end
    % retrieve index at angles 0 and pi
    [val1, ind1] = min( abs( angles ));
    [val2, ind2] = min( abs( angles - pi ));
    % Correlate images
    im1 = proj( :, rot_axis_roi2, ind1);
    im2 = flipud( proj( :, rot_axis_roi2, ind2) );
    PrintVerbose(verbose, '\n images to correlate: #%g at index %g and angle %g pi and #%g at index %g and angle %g pi', img_nums(ind1), ind1, val1 / pi, img_nums(ind2), ind2, (val2 + pi) / pi )
    out = ImageCorrelation( im1, im2, 0, 0);
    PrintVerbose(verbose, '\n respective shift: %g, %g', out.Xshift, out.Yshift)
    PrintVerbose(verbose, '\n calulated rotation axis offset: %g, %g', out.Xshift / 2)
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
        if vol_shape(1) <=  1
            vol_shape(1) = round( 1 + vol_shape(1) * (num_pix1 - 1) );
        end
        if vol_shape(2) <=  1
            vol_shape(2) = round( 1 + vol_shape(2) * (num_pix1 - 1) );
        end
        if vol_shape(3) <=  1
            vol_shape(3) = round( 1 + vol_shape(3) * (num_pix2 - 1) );
        end
    end
    if isempty( vol_size )
        vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];        
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
    subvol_shape = [vol_shape(1:2), num_sli];
    PrintVerbose(verbose, '\n rotation axis offset: %g', rotation_axis_offset );
    PrintVerbose(verbose, '\n rotation axis position: %g', num_pix1 / 2 + rotation_axis_offset );
    PrintVerbose(verbose, '\n GPU name: %s', gpu.Name );
    PrintVerbose(verbose, '\n GPU memory available: %g MiB (%.2f%%) of %g MiB', gpu.AvailableMemory/10^6, 100*gpu.AvailableMemory/gpu.TotalMemory, gpu.TotalMemory/10^6)
    PrintVerbose(verbose, '\n shape of reconstruction volume: [%g, %g, %g]', vol_shape )
    PrintVerbose(verbose, '\n shape of reconstruction subvolume: [%g, %g, %g]', subvol_shape )
    PrintVerbose(verbose, '\n memory required for reconstruction of a single slice: %g MiB (estimate)', rec_mem / 1024^2 )
    PrintVerbose(verbose, '\n memory of reconstructed volume: %g MiB', prod( vol_shape ) * 4 / 1024^2 )
    PrintVerbose(verbose, '\n maximum memory of subvolume: %g MiB', prod( subvol_shape ) * 4 / 1024^2 )
    PrintVerbose(verbose, '\n number of subvolume slabs: %g', num_slabs )
    PrintVerbose(verbose, '\n maximum number of slices per slab: %g', num_sli )    
    % Loop over slabs
    filt_direction = 1; % direction along the detector lines orthogonal to the rotation axis
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_phase_retrieval
        t = toc;   
        
        % Phase retieval filter
        edp = [energy, sample_detector_distance, bin*phys_pixel_size];
        [pf, pha_appendix] = PhaseFilter( phase_retrieval_method, binned_size, edp, reg_par, bin_filt, 'single');
          
        % reco phase dir
        if isempty( parfolder_reco )
            reco_phase_dir = [out_dir, filesep, 'reco_phase_', pha_appendix, filesep];
        else
            reco_phase_dir = [out_dir, filesep, 'reco_phase_', pha_appendix, filesep, parfolder_reco, filesep];
        end
        CheckAndMakePath( reco_phase_dir )
        PrintVerbose(verbose, '\n reco phase  : %s', reco_phase_dir)
        PrintVerbose(verbose, '\nPhase retrieval.')
        
                
        % Projections
        parfor nn = 1:num_proj
            pha = -real( ifft2( pf .* fft2( proj(:,:,nn) )) );
            
            proj(:,:,nn) = pha;
        end
        PrintVerbose(verbose, ' Elapsed time: %g s = %g min', toc-t, (toc-t)/60)        
end

%% Tomographic reco
if do_tomo
    PrintVerbose(verbose, '\n Tomographic reconstruction of')
    for nn = 1:num_slabs

        % Slice indices
        sli0 = 1 + ( nn - 1) * num_sli;
        sli1 = min( [nn * num_sli, vol_shape(3)] );

        % Filter sinogram
        PrintVerbose(verbose, ' \n  subvolume %g of %g: Filter sino: ', nn, num_slabs)                     
        filt = iradonDesignFilter('Ram-Lak', size(proj,1), 1);
        sino = real( ifft( bsxfun(@times, fft( NegLog(permute(proj(:,sli0:sli1,:), [1 3 2]), take_neg_log), [], 1), filt), [], 1, 'symmetric') );        
        PrintVerbose(verbose, 'done.')

        % Backprojection
        PrintVerbose(verbose, ' Backproject: ')    
        subvol_shape(3) = size( sino, 3);
        subvol_size = vol_size;
        subvol_size(5) = subvol_shape(3) / vol_shape(3) * vol_size(5);
        subvol_size(6) = subvol_shape(3) / vol_shape(3) * vol_size(6);
        vol = astra_parallel3D( sino, rot_global + angles, rotation_axis_offset, subvol_shape, subvol_size, pixel_size, link_data);
        PrintVerbose(verbose, 'done.')

        % Save subvolume
        if write_reco
            PrintVerbose(verbose, ' Save slices:')
            if do_phase_retrieval
                save_dir = reco_phase_dir;
            else
                save_dir = reco_dir;
            end
            CheckAndMakePath(reco_dir) 
            parfor ii = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_dir, sli0 + ii - 1);
                write32bitTIFfromSingle( filename, vol( :, :, ii) )
            end
            PrintVerbose(verbose, ' done.')
        end
    end
    PrintVerbose(verbose, '\n  Elapsed time: %g s = %g min', toc-t, (toc-t)/60 )

end

PrintVerbose(verbose, '\n Finished. Total time elapsed: %g s = %g min\n\n', toc, toc / 60 );
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
