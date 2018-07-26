function [proj, write] = p05_phase_retrieval( phase_retrieval, proj, write, verbose, visual_output)
% Phase retrieval for P05 data.
%
% Written by Julian Moosmann, last version: 2018-01-11
%
% [proj, reco_phase_path] = p05_phase_retrieval( phase_retrieval, edp, proj, write, verbose, visual_output)

t = toc;
PrintVerbose( verbose, '\nPhase retrieval: ')
PrintVerbose( verbose, '\n energy : %g eV', phase_retrieval.energy)
PrintVerbose( verbose, '\n sample detector distance : %g m', phase_retrieval.sample_detector_distance)
PrintVerbose( verbose, '\n pixel size : %g micron', phase_retrieval.eff_pixel_size_binned * 1e6)

method = phase_retrieval.method;
padding = phase_retrieval.padding;
reg_par = phase_retrieval.reg_par;
bin_filt = phase_retrieval.bin_filt;
cutoff_frequ = phase_retrieval.cutoff_frequ;

%% Phase retrieval filter
im_shape = [size(proj,1) , size(proj,2)];
im_shape_pad = (1 + padding) * im_shape;
edp = [phase_retrieval.energy, phase_retrieval.sample_detector_distance, phase_retrieval.eff_pixel_size_binned];
[phase_filter, pha_appendix] = PhaseFilter( method, im_shape_pad, edp, reg_par, bin_filt, cutoff_frequ, 'single');

% reco phase dir
if isempty( write.subfolder.reco )
    write.reco_phase_path = [write.path, filesep, 'reco_phase', filesep, pha_appendix, filesep];
else
    write.reco_phase_path = [write.path, filesep, 'reco_phase', filesep, pha_appendix, filesep, write.subfolder.reco, filesep];
end
CheckAndMakePath( write.reco_phase_path )
PrintVerbose( verbose, '\n reco_phase_path : %s', write.reco_phase_path)
PrintVerbose( verbose, '\n phase retrieval method : %s', method)

%% Retrieval
parfor nn = 1:size(proj, 3)
    % combined GPU and parfor usage requires memory management
    im = padarray( proj(:,:,nn), padding * im_shape, 'symmetric', 'post' );
    %im = padarray( gpuArray( proj(:,:,nn) ), raw_im_shape_binned, 'post', 'symmetric' );
    im = -real( ifft2( phase_filter .* fft2( im ) ) );
    proj(:,:,nn) = im(1:im_shape(1), 1:im_shape(2));
    %proj(:,:,nn) = gather( im(1:raw_im_shape_binned1, 1:raw_im_shape_binned2) );
end
pause(0.01)
PrintVerbose( verbose, '\n done in %g s (%.2f min)', toc-t, (toc-t)/60)

if visual_output
end

%% Post phase retrieval binning
phase_bin = phase_retrieval.post_binning_factor; % alias for readablity
if phase_bin > 1
    t = toc;
    PrintVerbose( verbose, '\nPost phase retrieval binning:')
    proj_bin = zeros( floor(size( proj ) ./ [phase_bin phase_bin 1] ), 'single');
    parfor nn = 1:size( proj, 3)
        proj_bin(:,:,nn) = Binning( proj(:,:,nn), phase_bin ) / phase_bin^2;
    end
    proj = proj_bin;
    PrintVerbose( verbose, ' done in %g s (%.2f min)', toc - t, (toc - t) / 60)
end

%% Save phase maps
if write.phase_map
    t = toc;
    PrintVerbose( verbose, '\nSave phase maps:')
    phase_map_path = [write.phase_map_path, pha_appendix, filesep];
    CheckAndMakePath( phase_map_path );
    % Save phase maps
    parfor nn = 1:size( proj, 3)
        filename = sprintf( '%sphase_%06u.tif', phase_map_path, nn);
        write32bitTIFfromSingle( filename, squeeze( rot90( proj( :, :, nn) ) ) )
    end
    pause(0.01)
    PrintVerbose( verbose, ' done in %.1f s (%.2f min)', toc - t, (toc - t) / 60)
end

%% Save sinograms
if write.phase_sino
    t = toc;
    PrintVerbose( verbose, '\nSave phase map sinogram:')
    sino_phase_path = write.sino_phase_path;
    CheckAndMakePath(sino_phase_path)
    % Save slices
    parfor nn = 1:size( proj, 2)
        filename = sprintf( '%ssino_%06u.tif', sino_phase_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    pause(0.01)
    PrintVerbose( verbose, ' done in %.1f s (%.2f min)', toc - t, (toc - t) / 60)
end
