function [proj, reco_phase_path] = p05_phase_retrieval( phase_retrieval_parameter, edp, proj, out_path, subfolder_reco, write, phase_map_path, verbose, visual_output)
% Phase retrieval for P05 data.
%
% Written by Julian Moosmann, last version: 2018-01-11
%
% [proj, reco_phase_path] = p05_phase_retrieval( phase_retrieval_parameter, edp, proj, subfolder_reco, write, phase_map_path, verbose, visual_output)

t = toc;
PrintVerbose( verbose, '\nPhase retrieval')
PrintVerbose( verbose, '\n energy : %g eV', edp(1))
PrintVerbose( verbose, '\n sample detector distance : %g m', edp(2))
PrintVerbose( verbose, '\n pixel size : %g micron', edp(3) * 1e6)

method = phase_retrieval_parameter.method;
padding = phase_retrieval_parameter.padding;
reg_par = phase_retrieval_parameter.reg_par;
bin_filt = phase_retrieval_parameter.bin_filt;
cutoff_frequ = phase_retrieval_parameter.cutoff_frequ;

%% Phase retrieval filter
im_shape = [size(proj,1) , size(proj,2)];
im_shape_pad = (1 + padding) * im_shape;
[phase_filter, pha_appendix] = PhaseFilter( method, im_shape_pad, edp, reg_par, bin_filt, cutoff_frequ, 'single');

% reco phase dir
if isempty( subfolder_reco )
    reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep];
else
    reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep, subfolder_reco, filesep];
end
CheckAndMakePath( reco_phase_path )
PrintVerbose( verbose, '\n reco_phase_path : %s', reco_phase_path)
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

%% Save phase maps
if write.phase_map
    t = toc;
    PrintVerbose( verbose, '\nSave phase maps:')
    phase_map_path = [phase_map_path, pha_appendix, filesep];
    CheckAndMakePath( phase_map_path );
    parfor nn = 1:size( proj, 3)
        filename = sprintf( '%sphase_%06u.tif', phase_map_path, nn);
        write32bitTIFfromSingle( filename, squeeze( rot90( proj( :, :, nn) ) ) )
    end
    pause(0.01)
    PrintVerbose( verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
end
