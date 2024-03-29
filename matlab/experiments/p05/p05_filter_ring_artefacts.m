function proj = pp_filter_ring_artefacts( ring_filter, proj, angles, par )
% Filter ring artefacts for P05 data.
%
% Written by Julian Moosmann.
%
% [proj, h] = pp_filter_ring_artefacts( ring_filter, proj, angles, par )

t = toc;
fprintf( '\nFilter ring artifacts using %s: ', ring_filter.method)
imsc1 = @(im) imsc( flipud( im' ) );
% sort angles for displaying and non-contiguous acquisition
[~, sorted_angle_index] = sort( angles );
sino_slice = round( size( proj, 2) / 2 );
sino_unfilt = squeeze( proj(:,sino_slice,sorted_angle_index) )';
switch lower( ring_filter.method )
    
    % Combined wavelet-FFT filter
    case 'wavelet-fft'
        
        dec_levels = ring_filter.waveletfft_dec_levels;
        wname = ring_filter.waveletfft_wname;
        sigma = ring_filter.waveletfft_sigma;
        
        parfor nn = 1:size( proj, 2)
            sino = squeeze( proj(:,nn,:) )';
            [d2, d1] = size( sino );
            sino = FilterStripesCombinedWaveletFFT( sino, dec_levels, wname, sigma )';
            sino = sino( 1:d1, 1:d2 );
            proj(:,nn,:) = sino;
        end
        sino_filt = squeeze( proj(:,sino_slice,sorted_angle_index) )';
        fprintf( ' done in %.1f s (%.2f min)', toc-t, (toc-t)/60)
        
        if par.visual_output
            figure('Name', 'Sinogram and ring filter', 'WindowState', 'maximized');
            
            subplot(3,1,1)
            imsc( sino_unfilt )
            axis equal tight
            title(sprintf('sino unfiltered, y = %u', sino_slice))
            colorbar
            
            subplot(3,1,2)
            imsc( sino_filt )
            axis equal tight
            title(sprintf('sino filtered, y = %u', sino_slice))
            colorbar
            
            subplot(3,1,3)
            imsc( sino_filt - sino_unfilt )
            axis equal tight
            title(sprintf('sino filt - sino, y = %u', sino_slice))
            colorbar
            
            drawnow
        end
        
    case 'jm'
        % Simple ring artifact filter
        if numel( ring_filter.jm_median_width ) > 1
            %% Combine if/else
            for nn = ring_filter.jm_median_width
                %m=sum(im);
                proj_mean = sum( proj, 3);
                
                %md=medfilt1(m,nn);
                proj_mean_med = medfilt2( proj_mean, [nn, 1], 'symmetric' );
                
                %f=md./m;
                mask = proj_mean_med ./ proj_mean;
                
                %im=bsxfun(@times,im,f);
                proj = bsxfun( @times, proj, mask);
            end
        else
            proj_mean = mean( proj, 3);
            proj_mean_med = medfilt2( proj_mean, [ring_filter.jm_median_width, 1], 'symmetric' );
            mask = proj_mean_med ./ proj_mean;
            proj = bsxfun( @times, proj, mask);
        end
        sino_filt = squeeze( proj(:,sino_slice,sorted_angle_index) )';
        
        fprintf( ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
        fprintf( '\n ring filter mask min/max: %f, %f', min( mask(:) ), max( mask(:) ) )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if par.visual_output
            figure('Name', 'Sinogram and ring filter', 'WindowState', 'maximized');
            
            subplot(2,2,1)
            imsc( sino_unfilt )
            axis equal tight
            title(sprintf('sino unfiltered, y = %u', sino_slice))
            colorbar %('Location', 'southoutside')
            
            subplot(2,2,2)
            imsc( sino_filt - sino_unfilt )
            axis equal tight
            title(sprintf('sino filt - sino, y = %u', sino_slice))
            colorbar %('Location', 'southoutside')
            
            subplot(2,2,3)
            imsc1( proj_mean )
            axis equal tight
            title(sprintf('mean projection'))
            colorbar %('Location', 'southoutside')
            
            subplot(2,2,4)
            imsc1( FilterHisto( mask, 5 ) )
            axis equal tight
            title(sprintf('mask for normalization'))
            colorbar %('Location', 'southoutside')
            
            drawnow
        end
end

