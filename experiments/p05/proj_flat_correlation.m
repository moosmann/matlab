function [proj, corr, toc_bytes] = proj_flat_correlation( proj, flat, image_correlation, par, roi_proj, roi_flat, toc_bytes )
% Correlate all projection with all flat-field to find the best matching
% pairs for the flat-field correction using method of choice.
%
% ARGUMENTS
% proj : stack of projections
% flat : stack of flat field
% image_correlation : struct, required fields see below
% par : struct, required fields see below
% roi_proj : projection ROI for correlation when offset shift is corrected
% roi_flat : reference ROI for correlation when offset shift is corrected
% 
% RETURNS
% proj : flat-field corrected projections
% corr : struct with correlation information
%
% Written by Julian Moosmann.
%
% [proj, corr] = proj_flat_correlation( proj, flat, image_correlation, par, roi_proj, roi_flat )


%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    roi_proj = [];
end
if nargin < 6
    roi_flat = [];
end

%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters in struct fields
method = image_correlation.method;
flat_corr_area1 = image_correlation.area_width;
flat_corr_area2 = image_correlation.area_height;
corr_num_flats = image_correlation.num_flats;
force_calc = image_correlation.force_calc;
im_shape_binned1 = image_correlation.im_shape_binned1;
im_shape_binned2 = image_correlation.im_shape_binned2;
flatcor_path = image_correlation.flatcor_path;
verbose = par.verbose;
x0 = par.offset_shift_x0;
%x1 = par.offset_shift_x1;
raw_bin = par.raw_bin;

im_shape_cropbin1 = size( proj, 1 );
num_proj_used = size( proj, 3);
num_ref_used = size( flat, 3);

%% CORRELATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
    
    case {'none', '', 'median'}
        fprintf( 'No correlation.' )
        corr = [];
        
    otherwise
        
        %% Load previously calculated correlation matrix & check if valid
        if ~force_calc
            % Loop over processed and scratch_cc
            filename_scra = sprintf( '%scorr.mat', flatcor_path);
            filename_proc = regexprep( filename_scra, 'scratch_cc', 'processed' );
            for filename = { filename_proc, filename_scra }
                % Load correlation matrix if exists
                if exist( filename{1} , 'file' )
                    load( filename{1},  'corr' );
                end
                % Check reco paramaters
                if exist( 'corr', 'var' )
                    force_calc = ~(...
                        isequal( corr.method, method) ...
                        && isequal( corr.size.proj, size( proj ) ) ...
                        && isequal( corr.size.flat, size( flat ) ) ...
                        && isequal( corr.raw_roi, image_correlation.raw_roi ) ...
                        && isequal( corr.raw_bin, image_correlation.raw_bin ) ...
                        && isequal( corr.proj_range, image_correlation.proj_range ) ...
                        && isequal( corr.ref_range, image_correlation.ref_range ) ...
                        );
                else
                    force_calc = 1;
                end
            end
        end
        
        % Correlate or not?
        if ~force_calc
            %% No correlation
            dt = '??';
            if isfield( corr, 'datetime' )
                dt = corr.datetime;
            end
            PrintVerbose( verbose, '\nLoading correlation computed on %s!', dt )
            
        else
            %% Correlate
            
            % Correlation ROI
            flat_corr_area1 = IndexParameterToRange( flat_corr_area1, im_shape_binned1 );
            flat_corr_area2 = IndexParameterToRange( flat_corr_area2, im_shape_binned2 );
            
            % ROI
            if isempty( roi_proj ) && isempty( roi_flat )
                roi_proj = proj(flat_corr_area1,flat_corr_area2,:);
                roi_flat = flat(flat_corr_area1,flat_corr_area2,:);
            end
            
            % Correlate flat fields: Compute correlation for each pair projection/flat-field
            PrintVerbose( verbose, '\nCorrelate projections & flat-fields. Method: %s.', method)
            t = toc;
            
            % Preallocation
            corr_mat = zeros( num_proj_used, num_ref_used );
            
            startS = ticBytes( gcp );
            switch lower( method )
                case 'entropy'
                    % entropy : entropy of ratio of p and f
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        %p = p(:);
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            f = p ./ f;
                            f = double( f );
                            corr_mat(pp,ff) = entropy( f );
                        end
                    end
                    
                case 'ssim-ml'
                    % ssim-ml: Matlab's structural similarity index (SSIM)
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            corr_mat(pp,ff) = - ssim( p, f );
                        end
                    end
                    
                case 'diff1-l1'
                    % anisotropic grad sum L1
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            d1 =  abs( abs( p ) - abs( f ) ) ;
                            corr_mat(pp,ff) = norm( d1(:), 1);
                        end
                    end
                    
                case 'diff1-l2'
                    % anisotropic grad sum L2
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            d1 =  abs( abs( p ) - abs( f ) ) ;
                            corr_mat(pp,ff) = norm( d1(:), 2);
                        end
                    end
                case 'diff2-l1'
                    % isotropic grad sum L1
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            d2 =  sqrt( abs( p.^2 - f.^2) );
                            corr_mat(pp,ff) = norm( d2(:), 1);
                        end
                    end
                case 'diff2-l2'
                    % isotropic grad sum L2
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            d2 =  sqrt( abs( p.^2 - f.^2) );
                            corr_mat(pp,ff) = norm( d2(:), 2);
                        end
                    end
                case 'std'
                    % std : standard deviation of ratio of p and f
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            f = p ./ f;
                            corr_mat(pp,ff) = std( f(:) );
                        end
                    end
                    
                case 'cross-entropy-12'
                    % cross entropy of proj and flat
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        p1 = imhist( normat( p(:) ) ) ./ numel( p );
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            p2 = imhist( normat( f(:) ), numel( p1 ) ) ./ numel( f );
                            m = boolean( (p1 == 0) + (p2 == 0) );
                            p1(m) = [];
                            p2(m) = [];
                            corr_mat(pp,ff) = - sum( p1 .* log( p2 ) );
                        end
                    end
                case 'cross-entropy-21'
                    % cross entropy of proj and flat
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        p1 = imhist( normat( p(:) ) ) ./ numel( p );
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            p2 = imhist( normat( f(:) ), numel( p1 ) ) ./ numel( f );
                            m = boolean( (p1 == 0) + (p2 == 0) );
                            p1(m) = [];
                            p2(m) = [];
                            corr_mat(pp,ff) = - sum( p2 .* log( p1 ) );
                        end
                    end
                    
                case 'cross-entropy-x'
                    % cross entropy of proj and flat
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        p1 = imhist( normat( p(:) ) ) ./ numel( p );
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            p2 = imhist( normat( f(:) ), numel( p1 ) ) ./ numel( f );
                            m = boolean( (p1 == 0) + (p2 == 0) );
                            p1(m) = [];
                            p2(m) = [];
                            corr_mat(pp,ff) = sum( p2 .* log2( p1 ) - p1 .* log2( p2 ) );
                        end
                    end
                    
                case 'cov'
                    % cov : cross covariance
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        p_mean = mean2( p );
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            f_mean = mean2( f );
                            cov_pf = mean2( ( p - p_mean ) .* (f - f_mean ) );
                            corr_mat(pp,ff) = - cov_pf;
                        end
                    end
                    
                case 'corr'
                    % corr : cross correlation
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        p_mean = mean2( p );
                        p_std = std( p(:) );
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            f_mean = mean2( f );
                            f_std = std( f(:) );
                            cov_pf = mean2( ( p - p_mean ) .* (f - f_mean ) );
                            corr_mat(pp,ff) = - cov_pf ./ ( p_std * f_std );
                        end
                    end
                    
                case 'ssim'
                    % ssim : own implementation of structural
                    % similarity index (SSIM) w/o Gaussian filter
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        p_mean = mean2( p );
                        p_std = std( p(:) );
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            f_mean = mean2( f );
                            f_std = std( f(:) );
                            cov_pf = mean2( ( p - p_mean ) .* (f - f_mean ) );
                            % Parameters
                            L = 1; % dynamic range of camera
                            c1 = ( 0.01 * L )^2;
                            c2 = ( 0.03 * L )^2;
                            c3 = c2 / 2;
                            % Components: luminance, contrast, structure
                            c_l = ( 2 * p_mean * f_mean + c1 ) / ( p_mean^2 + f_mean^2 + c1 );
                            c_c = ( 2 * p_std * f_std + c2 ) / ( p_std^2 + f_std^2 + c2 );
                            c_s = ( cov_pf + c3) / ( p_std * f_std + c3 );
                            corr_mat(pp,ff) = - c_l * c_c * c_s;
                        end
                    end
                    
                case 'ssim-g'
                    % ssim : own implementation of structural
                    % similarity index (SSIM) with Gaussian filter
                    radius = 1.5;
                    filtRadius = ceil( 3 * radius ); % 3 Standard deviations include >99% of the area.
                    filtSize = 2 * filtRadius + 1;
                    gaussFilterFcn = @(X) imgaussfilt( X, radius, 'FilterSize', filtSize, 'Padding','replicate');
                    parfor pp = 1:num_proj_used
                        p = roi_proj(:,:,pp);
                        p = gaussFilterFcn( p );
                        p_mean = mean2( p );
                        p_std = std( p(:) );
                        for ff = 1:num_ref_used
                            f = roi_flat(:,:,ff);
                            f = gaussFilterFcn( f );
                            f_mean = mean2( f );
                            f_std = std( f(:) );
                            cov_pf = mean2( ( p - p_mean ) .* (f - f_mean ) );
                            % Parameters
                            L = 1;% dynamic range of camera
                            c1 = ( 0.01 * L )^2;
                            c2 = ( 0.03 * L )^2;
                            c3 = c2 / 2;
                            % Components: luminance, contrast, structure
                            c_l = ( 2 * p_mean * f_mean + c1 ) / ( p_mean^2 + f_mean^2 + c1 );
                            c_c = ( 2 * p_std * f_std + c2 ) / ( p_std^2 + f_std^2 + c2 );
                            c_s = ( cov_pf + c3) / ( p_std * f_std + c3 );
                            corr_mat(pp,ff) = - c_l * c_c * c_s;
                        end
                    end
                    
            end % switch lower( method )
            toc_bytes.correlation = tocBytes( gcp, startS );
            
            % Save correlation matrix
            corr.mat = corr_mat;
            corr.method = method;
            corr.size.proj = size( proj );
            corr.size.flat = size( flat );
            corr.datetime = datetime;
            corr.raw_roi = image_correlation.raw_roi;
            corr.raw_bin = image_correlation.raw_bin;
            corr.proj_range = image_correlation.proj_range;
            corr.ref_range = image_correlation.ref_range;
            CheckAndMakePath( flatcor_path )
            save( sprintf( '%scorr.mat', flatcor_path), 'corr' )
            
            PrintVerbose( verbose, ' Done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
            
        end %  if ~force_calc
        
end % switch method

%% Plot correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.visual_output
    figure( 'Name', 'Correlation of projections and flat-fields');
    mid = round( num_proj_used / 2 );
    f = @(mat,pp) mat(pp,:)';
    Y = [ f( corr.mat, 1), f( corr.mat, mid), f( corr.mat, num_proj_used) ];
    plot( Y, '.' )
    legend( 'first proj', 'mid proj', 'last proj' )
    axis tight
    title(sprintf('correlation method: %s', method))
    xlabel( 'flat field index' )
    ylabel( 'measure' )
    drawnow
end

%% CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
    case {'none', '', 'median'}
        
        % Flat field correction without correlation
        PrintVerbose(verbose, '\nFlat-field correction w/o correlation.')
        
        flat_median = median( flat, 3);
        parfor nn = 1:num_proj_used
            im = proj(:, :, nn);
            flat_median_shifted = circshift( flat_median, round( -x0(nn) / raw_bin ), 1 );
            proj(:, :, nn) = im ./ flat_median_shifted;
        end
        
    otherwise
        PrintVerbose(verbose, '\nFlat-field correction with correlation.')
        t = toc;
        
        %         % Memory requirement and parpool processing
        %         [mem_free, ~, mem_total] = free_memory;
        %         mem_req_tot = Bytes( proj ) + Bytes(flat) * (poolsize + 1);
        %         mem_req_par = Bytes(flat) * poolsize;
        %         if mem_req_tot < 0.75 * mem_total && mem_req_par < mem_free
        %             pflag = 1;
        %             fprintf( ' Use parpool.' )
        %         else
        %             pflag = 0;
        %             fprintf( ' No parpool.' )
        %        end
        
        [~, corr_mat_pos] = sort( normat( corr.mat ), 2);
        % Flat field correction
        flat_ind = corr_mat_pos(:,1:corr_num_flats);
        
        %         if pflag
        %
        %             %% prallel
        %             startS = ticBytes( gcp );
        
        for nn = 1:num_proj_used
            
            % Mean over flat fields
            flat_ind_nn = flat_ind(nn,:);
            flat_mean = mean( flat(:,:,flat_ind_nn), 3);
            
            % Binned shift (shift, not first pixel)
            shift = ( x0(nn) - 1 ) / raw_bin;
            shift_int = floor( shift );
            shift_sub = shift - shift_int;
            
            % shift flat
            if mod( shift_sub, 1 ) ~= 0
                % crop flat at integer shift, then shift subpixel
                xx = shift_int + (1:im_shape_cropbin1+1);
                flat_mean_shifted = imtranslate( flat_mean(xx,:), [0 -shift_sub], 'linear' );
            else
                xx = shift_int + (1:im_shape_cropbin1);
                flat_mean_shifted = flat_mean(xx,:);
            end
            
            % flat field correction
            p = proj(:, :, nn);
            p = p ./ flat_mean_shifted(1:im_shape_cropbin1,:) ;
            proj(:,:,nn) = p;
            %imsc( p(1:200,:) ),axis equal tight
            
        end
        %            toc_bytes.correction = tocBytes( gcp, startS );
        
        PrintVerbose( verbose, ' Done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
end
