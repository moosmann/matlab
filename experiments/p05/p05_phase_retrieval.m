function [proj, write, tint] = p05_phase_retrieval( proj, phase_retrieval, tomo, write, interactive_mode, verbose, visual_output)
% Phase retrieval for P05 data.
%
% Written by Julian Moosmann, 2018-01-11, last version: 2018-07-27
%
% [proj, write] = p05_phase_retrieval( proj, phase_retrieval, tomo, write, interactive_mode, verbose, visual_output)

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = toc;
PrintVerbose( verbose, '\nPhase retrieval: ')

method = phase_retrieval.method;
padding = phase_retrieval.padding;
reg_par = phase_retrieval.reg_par;
bin_filt = phase_retrieval.bin_filt;
cutoff_frequ = phase_retrieval.cutoff_frequ;
edp = [phase_retrieval.energy, phase_retrieval.sample_detector_distance, phase_retrieval.eff_pixel_size_binned];
if edp(1) <= 0
    error( 'Energy is zero!' );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interactive mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tint = 0;
if exist( 'interactive_mode', 'var' ) && isfield( interactive_mode, 'phase_retrieval' ) && interactive_mode.phase_retrieval
    tint = toc;
    plot_ctf_mask(1) = 0;
    
    if tomo.rot_axis.tilt ~= 0
        error( 'Rotation axis tilt not supported in interactive phase retrieval mode.' )
    end
    
    if plot_ctf_mask
        proj1 = proj(:,:,1);
        proj_mean = mean( proj, 3);
    end
            
    fprintf( '\n\nENTER INTERACTIVE PHASE RETRIEVAL MODE' )
    
    reg_par_def_range = 0:6;
    
    if interactive_mode.slice_number > 1
        slice = interactive_mode.slice_number;
    elseif interactive_mode.slice_number <= 1 && interactive_mode.slice_number >= 0
        slice = round( ( size( proj, 2) - 1) * interactive_mode.slice_number + 1 );
    else
        slice = ceil( size( proj, 2) / 2 );
    end
    
    % Tomography parameters
    vol_shape = assign_from_struct( tomo, 'vol_shape', [] );
    vol_size = assign_from_struct( tomo, 'vol_size', [] );
    butterworth_filtering = assign_from_struct( tomo.butterworth_filter, 'apply', 0 );
    tomo.rot_axis.offset = tomo.rot_axis.offset + tomo.offset_shift + eps;
    [num_pix, num_row, num_proj] = size( proj );
       
    %  Size, shape
    [vol_shape, vol_size] = volshape_volsize( proj, vol_shape, vol_size, tomo.rot_axis.offset, 0);
    if isempty( vol_shape )
        vol_shape = [num_pix, num_pix, 1];
    else
        vol_shape(3) = 1;
    end
    tomo.vol_shape = vol_shape;
    if isempty( vol_size )
        vol_size = [-num_pix/2 num_pix/2 -num_pix/2 num_pix/2 -0.5 0.5];
    else
        vol_size(5) = -0.5;
        vol_size(6) = 0.5;
    end
    tomo.vol_size = vol_size;
    
    if strcmpi( tomo.algorithm, 'fbp' )
        % Ramp filter
        fbp_filt = iradonDesignFilter('Ram-Lak', 2 * num_pix, 1);
        
        % Butterworth filter
        if butterworth_filtering
            [b, a] = butter(1, 0.5);
            bw = freqz(b, a, numel(fbp_filt) );
            fbp_filt = fbp_filt .* bw;
        end
    end
    
    % Metrics
    mask_rad = 0.95;
    mask_val = 0;
    reco_metric(1).name = 'mean';
    reco_metric(2).name = '1-norm';    
    reco_metric(3).name = 'neg';
    reco_metric(4).name = 'iso-grad';
    reco_metric(5).name = 'laplacian';
    reco_metric(6).name = 'entropy';
    reco_metric(7).name = 'entropyML';    
    reco_metric(8).name = '2-norm';
    
    %% Interaction loop
    first_round = 1;
    slice_old = 0;
    reg_par = reg_par_def_range;
    slab_size = 100;
    while ~isscalar( reg_par ) || first_round
        
        % Print parameters
        fprintf( '\n current reg_par range : %g', reg_par(1) )
        fprintf( ',%g', reg_par(2:end) )
        fprintf( '\n current slice : %u', slice )
        fprintf( '\n current method : %s', method )
        fprintf( '\n current [energy, distance, pixelsize] : [%.1f keV, %.1f cm, %.3f mu]', edp(1)/1000, edp(2)*1e2, edp(3) * 1e6 )
        if sum( strcmp( method, {'qpcut'} ) )
            fprintf( '\n current cut-off frequency : %.2f', cutoff_frequ )
        end
        if sum( strcmp( method, {'qp', 'qp2', 'qpcut'} ) )
            fprintf( '\n current binary filter threshold : %.3f', bin_filt )
        end
        fprintf( '\n current padding factor : %u', padding )
        fprintf( '\n current vertical slab size : %u', slab_size )
        
        if plot_ctf_mask
            % FFT
            epsilon = 0;
            fun = @(im) Binning( FilterHisto(fftshift(log( 10^(-epsilon) + abs( fft2( padarray( im , padding*size(im ), 'symmetric', 'post' ) ) ) ) ), 3), 2*padding);
            proj1_fft = fun( proj1 );
            proj_mean_fft = fun( proj_mean );
            
            % CTF binary mask
            [~, mask_half] = BinaryMaskCTF( size(proj1_fft), edp(1), edp(2), ( 1 + 0*padding ) * edp(3), 1, 0.1);
            
            % Image plus half mask
            proj1_fft_mask = normat( proj1_fft ) .* mask_half ;
            proj_mean_fft_mask = normat( proj_mean_fft ) .* mask_half;
            
            figure( 'Name', 'CHECK PHASE RETRIEVAL PARAMETERS: compare Fourier transformed projection and CTF mask' )
            
            subplot(1,2,1)
            imsc( proj1_fft_mask )
            axis equal tight
            title(sprintf('FFT projection #1 * CTF mask'))
            
            subplot(1,2,2)
            imsc( proj_mean_fft_mask )
            axis equal tight
            title(sprintf('FFT mean projection * CTF mask'))
            
            drawnow
        end
        
        if ~first_round
            
            % Metric preallocation
            for nn = 1:numel( reco_metric )
                reco_metric(nn).val = zeros( numel(reg_par), 1);
            end
            
            %% Reconstruction loop overregularization parameter
            pha = zeros( [vol_shape(1), vol_shape(2), numel( reg_par )], 'single' );
            for nn = 1:numel( reg_par )
                
                fprintf( '\n  Step %2u: ', nn )
                
                % Slab, pad, FT                
                dz = round( slab_size / 2 );
                rows = slice + (-dz:dz);
                rows = rows - max( 0, max( rows(:) - size( proj, 2) ) );
                num_row_slab = numel( rows );
                im_shape = [num_pix, num_row_slab, num_proj];
                if ~isequal( slice, slice_old )
                    t = toc;
                    fprintf( ' fft o pad: ' )
                    slab_ft = fft2( padarray( proj(:, rows, :), [padding,padding,0] .* im_shape, 'symmetric', 'post' ) );
                    fprintf( ' %.1f s.', toc - t)
                end
                
                % Phase retrieval filter
                im_shape_pad = (1 + padding) * im_shape;
                phase_filter = PhaseFilter( method, im_shape_pad, edp, reg_par(nn), bin_filt, cutoff_frequ, 'single');
                
                % Retrieval plus FBP filter
                t = toc;
                fprintf( ' filter:' )
                tmp = bsxfun( @times, slab_ft, phase_filter );
                if strcmpi( tomo.algorithm, 'fbp' )
                    tmp = bsxfun( @times, tmp, fbp_filt );
                end
                fprintf( ' %.1f s.', toc - t)
                t = toc;
                fprintf( ' real o ifft:' )
                slab = -real( ifft2( tmp ) );
                sino = slab(1:num_pix,dz + 1,:);
                fprintf( ' %.1f s.', toc - t)
                
                % Tomography
                t = toc;
                fprintf( ' tomo:' )
                reco = astra_parallel2D( tomo, permute( sino, [3 1 2]) );
                %im = FilterOutlier( im );
                pha(:,:,nn) = reco;
                fprintf( ' %.1f s.', toc - t)
                
                % Metrics
                t = toc;
                fprintf( ' metrics:' )
                reco = double( MaskingDisc( reco, mask_rad, mask_val) ) * 2^16;
                % mean
                reco_metric(1).val(nn) = mean2( reco );
                % mean abs
                reco_metric(2).val(nn) = mean2( abs( reco ) );
                % mean negativity
                reco_metric(3).val(nn) = - mean( reco( reco <= 0 ) );
                % isotropic gradient
                [g1, g2] = gradient(reco);
                reco_metric(4).val(nn) = mean2( sqrt( g1.^2 + g2.^2 ) );
                % laplacian
                reco_metric(5).val(nn) = mean2( abs( del2( reco ) ) );
                % entropy
                p = histcounts( reco(:) );
                p = p(p>0);
                p = p / sum( p );
                reco_metric(6).val(nn) = -sum( p .* log2( p ) );
                % entropy built-in
                reco_metric(7).val(nn) = -entropy( reco );
                % 2-norm
                reco_metric(8).val(nn) = norm( reco(:), 2 );
                fprintf( ' %.1f s.', toc - t)

            end
            
            % Normalize for ease of plotting and comparison
            for nn = 1:numel(reco_metric)
                reco_metric(nn).val = normat( reco_metric(nn).val );
            end
            
            % Metric minima
            [~, min_pos] = min(cell2mat({reco_metric(:).val}));
            [~, max_pos] = max(cell2mat({reco_metric(:).val}));
            
            % Print image number, regularization parameter, and different metrics
            fprintf( '\n')
            fprintf( '%10s', 'image no.', 'reg par', 'reg val', reco_metric.name)
            for nn = 1:numel(reg_par)
                if reg_par(nn) == phase_retrieval.reg_par
                    cprintf( 'Green', sprintf('\n%10u%10.2f%10.2g', nn, reg_par(nn), 10^-reg_par(nn) ) )
                else
                    cprintf( 'Black', '\n%10u%10.2f%10.2g', nn, reg_par(nn), 10^-reg_par(nn) )
                end
                
                for mm = 1:numel(reco_metric)
                    if min_pos(mm) == nn
                        cprintf( 'Red', '%10.3g', reco_metric(mm).val(nn) )
                    elseif max_pos(mm) == nn
                        cprintf( 'Blue', '%10.3g', reco_metric(mm).val(nn) )
                    else
                        cprintf( 'Black', '%10.3g', reco_metric(mm).val(nn) )
                    end
                end
            end
            
            % Plot metrics
            figure('Name', 'REGULARIZATION PARAMETER: metrics');
            x = 1:numel( reco_metric );
            Y = cell2mat({reco_metric(x).val});
            plot( reg_par, Y, '-+');            
            axis tight
            xlabel( 'regularization parameter' )
            legend( reco_metric(x).name )
            ax1 = gca;            
            ax2 = axes( 'Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
            line(1:numel( reg_par ), 0, 'Parent', ax2 )
            xlabel( 'index (image no.)' )
            set( ax1, 'YTick', [] ) % 'XTickMode', 'auto', 'XMinorTick', 'on')
            set( ax2, 'YTick', [] )
            title(sprintf('rotation axis: metrics VS regularization parameter'))
            drawnow
            
            % Play
            nimplay( pha, 1, [], 'PHASE RETRIEVAL: sequence of reconstructed slices using different phase retrieval parameters')
            
        else
            first_round = 0;
        end
        
        % Input
        reg_par = '';
        while ischar( reg_par )
            reg_par = input( '\n\nENTER REG_PAR OR A RANGE OF REG_PAR (if empty: use default range, scalar: use value & exit loop, ''s'': change slice number,\n ''m'':change method, ''b'': binary filter threshold, ''c'': change cut-off frequency, ''p'': padding, ''d'': debug mode): ');
            if isempty( reg_par )
                reg_par = reg_par_def_range;
            else
                if ischar( reg_par )
                    switch reg_par
                        case 's'
                            slice_old = slice;
                            slice = input( sprintf( '\n\nENTER ABSOLUTE [1,%u] OR RELATIVE [0,1] SLICE NUMBER : ', size(proj,2)) );
                            if slice <= 1 && slice >= 0
                                slice = round( (num_row - 1) * slice + 1 );
                            end
                            fprintf( ' \n new slice : %u', slice );
                        case 'edp'
                            edp = input( '\n\nENTER [energy/keV, distance/m, pixelsize/m] : ' );
                        case 'm'
                            method = input( '\n\nENTER METHOD (''tie'',''qp'',''qp2'',''ctf'',''qpcut'') : ' );
                        case 'c'
                            cutoff_frequ = input( '\n\nENTER CUT-OFF FREQUENCY (in rad) : ' );
                        case 'b'
                            bin_filt = input( '\n\nENTER BINARY FILTER THRESHOLD (in [0,1], typically 0.01-0.1) : ' );
                        case 'p'
                            padding = input( '\n\nENTER PADDING FACTOR (integer 0,1,2,..) : ' );
                        case 'ss'
                            slab_size = input( '\n\nENTER SLAB SIZE (integer>1, typically 50-200) : ' );
                        case 'd'
                            keyboard
                    end % switch
                end % ischar
            end % if ~isempty( reg_par )
        end % of while loop: ischar( reg_par )
        
    end % of while loop ~isscalar( reg_par )
    
    phase_retrieval.method = method;
    phase_retrieval.reg_par = reg_par;
    phase_retrieval.cutoff_frequ = cutoff_frequ;
    phase_retrieval.bin_filt = bin_filt;
    phase_retrieval.padding = padding;
    phase_retrieval.energy = edp(1);
    phase_retrieval.sample_detector_distance = edp(2);
    phase_retrieval.eff_pixel_size_binned = edp(3);
    
    tint = toc - tint;
    fprintf( '\nEND OF INTERACTIVE MODE' )  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Phase retrieval filter
im_shape = [size(proj,1) , size(proj,2)];
im_shape_pad = (1 + padding) * im_shape;
[phase_filter, pha_appendix] = PhaseFilter( method, im_shape_pad, edp, reg_par, bin_filt, cutoff_frequ, 'single');

% reco phase dir
if isempty( write.subfolder.reco )
    write.reco_phase_path = [write.path, filesep, 'reco_phase', filesep, pha_appendix, filesep];
else
    write.reco_phase_path = [write.path, filesep, 'reco_phase', filesep, pha_appendix, filesep, write.subfolder.reco, filesep];
end
CheckAndMakePath( write.reco_phase_path )
PrintVerbose( verbose, '\n energy : %g eV', phase_retrieval.energy)
PrintVerbose( verbose, '\n sample detector distance : %g m', phase_retrieval.sample_detector_distance)
PrintVerbose( verbose, '\n pixel size : %g micron', phase_retrieval.eff_pixel_size_binned * 1e6)
PrintVerbose( verbose, '\n phase retrieval method : %s', method)
PrintVerbose( verbose, '\n reco_phase_path : %s', write.reco_phase_path)

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
PrintVerbose( verbose, '\n done in %g s (%.2f min)', toc-t-tint, (toc-t-tint)/60)

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
