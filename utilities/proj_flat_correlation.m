function [proj, corr] = proj_flat_correlation( proj, flat, correlation_method, flat_corr_area1, flat_corr_area2, raw_im_shape_binned, corr_shift_max_pixelshift, corr_num_flats, decimal_round_precision, flatcor_path, verbose, visualOutput)
% Correlate all projection with all flat-field to find the best matching
% pairs for the flat-field correction using method of choice.
%
% correlation_method :
%   'shift' : cross correlation equals convolution of complex conjugate of
%   p(-x) i.e. rot180(p) and f(x), (p^*)(-x) ** f(x)
%
%
% Written by Julian Moosmann. First version: 2017-08-23. Last modifcation:
% 2017-08-23

%% Projection / flat-field correlcation
switch correlation_method

    case {'none', ''}
        % Flat field correction without correlation
        flat_median = median( flat, 3);
        proj = bsxfun( @times, proj, flat_median);
        corr = [];        
        
    otherwise
        
        num_proj_used = size( proj, 3);
        num_ref_used = size( flat, 3);
        
        % Correlate flat fields
        PrintVerbose(verbose, '\nCorrelate projections and flat-fields.')
        t = toc;
        
        % Correlation ROI
        flat_corr_area1 = IndexParameterToRange(flat_corr_area1, raw_im_shape_binned(1));
        flat_corr_area2 = IndexParameterToRange(flat_corr_area2, raw_im_shape_binned(2));
        flat_roi = flat(flat_corr_area1, flat_corr_area2, :);
        proj_roi = proj(flat_corr_area1, flat_corr_area2, :);
        
        % Preallocation
        foo = zeros( num_proj_used, num_ref_used);
        c_shift_1 = foo;
        c_shift_2 = foo;
        c_diff1_l1 = foo;
        c_diff1_l2 = foo;
        c_diff2_l1 = foo;
        c_diff2_l2 = foo;
        c_std = foo;
        c_ent = foo;
        c_cross_entropy12 = foo;
        c_cross_entropy21 = foo;
        c_cross_entropyx = foo;
        c_cov = foo;
        c_corr = foo;
        c_ssim = foo;
        c_ssim_ml = foo;
        c_l = foo;
        c_c = foo;
        c_s = foo;
        
        % Compute shift/correlation for each pair projection/flat-field
        for ff = 1:num_ref_used
            
            % flat field
            f = flat_roi(:,:,ff);
            
            parfor pp = 1:num_proj_used
                
                % projection
                p = proj_roi(:,:,pp);
                
                switch correlation_method
                    
                    case 'shift'
                        % shift via image cross correlation
                        out = ImageCorrelation( p, f, 0, 0, 0, 0, 1);
                        c_shift_1(pp,ff) = round( out.shift1, decimal_round_precision );
                        c_shift_2(pp,ff) = round( out.shift2, decimal_round_precision) ; % relevant shift
                        
                    case 'diff'
                        
                        % differences
                        d1 =  abs( abs( p ) - abs( f ) ) ;
                        d2 =  sqrt( abs( p.^2 - f.^2) );
                        
                        % anisotropic grad sum L1
                        c_diff1_l1(pp,ff) = norm( d1(:), 1);
                        
                        % anisotropic grad sum L2
                        c_diff1_l2(pp,ff) = norm( d1(:), 2);
                        
                        % isotropic grad sum L1
                        c_diff2_l1(pp,ff) = norm( d2(:), 1);
                        
                        % isotropic grad sum L2
                        c_diff2_l2(pp,ff) = norm( d2(:), 2);
                        
                        
                    case 'std'
                        % std : standard deviation of ratio of p and f
                        c_std(pp,ff) = std2( p ./ f );
                        
                    case 'entropy'
                        % entropy : entropy of ratio of p and f
                        c_ent(pp,ff) = entropy( double( p ./ f ) );
                        
                    case {'cross-entropy12', 'cross-entropy21', 'cross-entropyx'}
                        % cross entropy of proj and flat
                        p1 = imhist( normat( p(:) ) ) ./ numel( p );
                        p2 = imhist( normat( f(:) ), numel( p1 ) ) ./ numel( f );
                        m = boolean( (p1 == 0) + (p2 == 0) );
                        p1(m) = [];
                        p2(m) = [];
                        c_cross_entropy12(pp,ff) = - sum( p1 .* log( p2 ) );
                        c_cross_entropy21(pp,ff) = - sum( p2 .* log( p1 ) );
                        c_cross_entropyx(pp,ff) = sum( p2 .* log2( p1 ) - p1 .* log2( p2 ) );
                        
                    case {'cov', 'corr', 'ssim'}
                        
                        % Inputs
                        p_mean = mean2( p );
                        p_std = std2( p );
                        f_mean = mean2( f );
                        f_std = std2( f );
                        cov_pf = mean2( ( p - p_mean ) .* (f - f_mean ) );
                        
                        % cov : cross covariance
                        c_cov(pp,ff) = - cov_pf;
                        
                        % corr : cross correlation
                        c_corr(pp,ff) = - cov_pf ./ ( p_std * f_std );
                        
                        
                        % ssim : structural similarity index (SSIM)
                        
                        % Dynamic range of camera
                        %                     switch lower( cam )
                        %                         case 'kit'
                        %                             L = 2^12;
                        %                         case 'ehd'
                        %                             L = 2^16;
                        %                     end
                        % L = round(max(max(flat_roi(:)),max(flat_roi(:))) - min(min(flat_roi(:)),min(flat_roi(:))));
                        L = 1;
                        
                        % Parameters
                        c1 = ( 0.01 * L )^2;
                        c2 = ( 0.03 * L )^2;
                        c3 = c2 / 2;
                        
                        % Components: luminance, contrast, structure
                        c_l(pp,ff) = ( 2 * p_mean * f_mean + c1 ) / ( p_mean^2 + f_mean^2 + c1 );
                        c_c(pp,ff) = ( 2 * p_std * f_std + c2 ) / ( p_std^2 + f_std^2 + c2 );
                        c_s(pp,ff) = ( cov_pf + c3) / ( p_std * f_std + c3 );
                        
                        c_ssim(pp,ff) = - c_l(pp,ff) * c_c(pp,ff) * c_s(pp,ff);
                        
                    case 'ssim-ml'
                        % Matlab's structural similarity index (SSIM)
                        c_ssim_ml(pp,ff) = - ssim( proj_roi(:,:,pp), f ); %'DynamicRange', 'uint16'
                end
            end
        end
        
        switch correlation_method
            case 'shift'
                corr_mat = c_shift_2;
            case 'diff'
                corr_mat = c_diff1_l2;                
            case 'std'
                corr_mat = c_std;
            case 'entropy'
                corr_mat = c_ent;
            case 'cross-entropy12'
                corr_mat = c_cross_entropy12;
            case 'cross-entropy21'
                corr_mat = c_cross_entropy21;
            case 'cross-entropyx'
                corr_mat = c_cross_entropyx;
            case 'ssim'
                corr_mat = c_ssim;
            case 'ssim-ml'
                corr_mat = c_ssim_ml;
            case 'cov'
                corr_mat = c_cov;
            case 'corr'
                corr_mat = c_corr;
        end
        
        if strcmp( correlation_method, 'all' )
            % sorted measures: position and values
            [corr.diff1_l1_val,corr.diff1_l1_pos] = sort( c_diff1_l1, 2 );
            [corr.diff1_l2_val,corr.diff1_l2_pos] = sort( c_diff1_l2, 2 );
            [corr.diff2_l1_val,corr.diff2_l1_pos] = sort( c_diff2_l1, 2 );
            [corr.diff2_l2_val,corr.diff2_l2_pos] = sort( c_diff2_l2, 2 );
            [corr.std_val,corr.std_pos] = sort( c_std, 2 );
            [corr.ent_val,corr.ent_pos] = sort( c_ent, 2 );
            [corr.cross_ent12_val,corr.cross_ent12_pos] = sort( c_cross_entropy12, 2 );
            [corr.cross_ent21_val,corr.cross_ent21_pos] = sort( c_cross_entropy21, 2 );
            [corr.cross_entx_val,corr.cross_entx_pos] = sort( c_cross_entropyx, 2 );
            [corr.cov_val,corr.cov_pos] = sort( c_cov, 2 );
            [corr.corr_val,corr.corr_pos] = sort( c_corr, 2 );
            [corr.ssim_val,corr.ssim_pos] = sort( c_ssim, 2 );
            [corr.ssim_ml_val,corr.ssim_ml_pos] = sort( c_ssim_ml, 2 );
        end
        
        PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        
        %% Visual output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if visualOutput(1)
            figure('Name', 'Correlation of projections and flat-fields');
            
            switch correlation_method
                case 'shift'
                    [~, flat_corr_shift_min_pos_x] =  min ( abs( c_shift_1), [], 2);
                    [~, flat_corr_shift_min_pos_y] =  min ( abs( c_shift_2), [], 2);
                    
                    m = 2; n = 1;        
                    
                    subplot(m,n,1);
                    Y = abs(arrayfun(@(x) (c_shift_2(x,flat_corr_shift_min_pos_y(x))), 1:num_proj_used));
                    plot(Y, '.')
                    axis  tight
                    title(sprintf('minimal absolute vertical shift along rotation axis'))
                    
                    subplot(m,n,2);
                    plot( arrayfun(@(x) (c_shift_1(x,flat_corr_shift_min_pos_x(x))), 1:num_proj_used) ,'.')
                    axis  tight
                    title(sprintf('minimal absolute horizontal shift (unused)'))
                    
                case {'cov', 'corr', 'ssim'}
                    
                    m = 2; n = 1;
                    subplot(m,n,1);                    
                    f = @(x) normat(x(1,:))';
                    Y = [ f(c_cov), f(c_corr), f(c_ssim)];
                    plot( Y, '-' )
                    legend( 'cov', 'corr', 'ssim' )
                    axis tight
                    title(sprintf('correlation measures: projection #1'))
                    
                    subplot(m,n,2);
                    f = @(x) normat(x(1,:))';
                    Y = [ f(c_ssim), f(-c_l), f(-c_c), f(-c_s) ];
                    plot( Y, '-' )
                    legend( 'ssim', 'luminance', 'contrast', 'structure' )
                    axis tight
                    title(sprintf('SSIM and components: projection #1'))
                
                case 'all'
                    
%                     m = 2; n = 1;
%                     subplot(m,n,1);
%                     %f = @(x) normat( min( x, [], 2));
%                     f = @(x) normat(x(1,:))';
%                     Y = [f(c_diff1_l2), f(c_std), f(c_ent), f(c_cross_entropy12), f(c_cross_entropy21), f(c_cross_entropyx), f(c_cov), f(c_corr), f(c_ssim), f(c_ssim_ml)];
%                     plot( Y, '-' )
%                     legend('iso diff L2', 'ratio std dev', 'ratio entropy', 'cross entropy 12', 'cross entropy 21','cross entropy x', 'cov', 'corr', 'ssim', 'ssim-ml' )
%                     axis tight
%                     title(sprintf('correlation measures: projection #1'))

                    
                otherwise
                                        
                    mid = round( num_proj_used / 2 );
                    f = @(mat,pp) mat(pp,:)';
                    Y = [ f( corr_mat, 1), f( corr_mat, mid), f( corr_mat, num_proj_used) ];
                    plot( Y, '.' )
                    legend( 'first proj', 'mid proj', 'last proj' )
                    axis tight
                    title(sprintf('correlation method: %s', correlation_method))      
                    xlabel( 'flat field index' )
                    ylabel( 'measure' )
                    
            end
            drawnow
        end
        
        
        %% Flat- and dark-field correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        PrintVerbose(verbose, '\nFlat- and dark-field correction.')
        t = toc;
        switch correlation_method
            case 'shift'
                % best match
                if corr_shift_max_pixelshift == 0                    
                    [~, pos] = min( abs( corr_mat ), [], 2 );
                    parfor nn = 1:num_proj_used
                        proj(:, :, nn) = proj(:, :, nn) ./ flat(:, :, pos(nn));
                    end
                    
                    % use all flats which are shifted less pixels than corr_shift_max_pixelshift
                elseif corr_shift_max_pixelshift > 0
                    nflats = zeros(1, num_proj_used);
                    parfor nn = 1:num_proj_used
                        vec = 1:num_ref_used;
                        flat_ind = vec( abs( c_shift_2(nn, :) ) < corr_shift_max_pixelshift );
                        if numel( flat_ind ) > corr_num_flats
                            flat_ind( corr_num_flats + 1:end ) = [];
                        end
                        if isempty( flat_ind )
                            [~, flat_ind] = min( abs( c_shift_2(nn, :) ) );
                        end
                        nflats(nn) = numel(flat_ind);
                        proj(:, :, nn) = proj(:, :, nn) ./ squeeze( mean( flat(:, :, flat_ind), 3) );
                    end
                    PrintVerbose(verbose, '\n number of flats used per projection: [mean, min, max] = [%g, %g, %g]', mean( nflats ), min( nflats ), max( nflats) )
                else
                    error('Value of maximum shift (%g) is not >= 0', corr_shift_max_pixelshift)
                end
                
            otherwise
                
                [corr_mat_val, corr_mat_pos] = sort( normat( corr_mat ), 2);
                corr.mat = corr_mat;
                corr.mat_val = corr_mat_val;
                corr.mat_pos = corr_mat_pos;
                
                % Save correlation matrix
                CheckAndMakePath( flatcor_path )
                save( sprintf( '%s/corr_mat_val.mat', flatcor_path), 'corr_mat_val' )
                save( sprintf( '%s/corr_mat_pos.mat', flatcor_path), 'corr_mat_pos' )
                
                % Flat field correction
                flat_ind = corr_mat_pos(:,1:corr_num_flats);
                parfor nn = 1:num_proj_used
                    proj(:, :, nn) = proj(:, :, nn) ./ squeeze( mean( flat(:, :, flat_ind(nn,:)), 3) );
                end
                
        end
        PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
end
