function [proj, corr, toc_bytes] = proj_flat_correlation( proj, flat, image_correlation, par, write, roi_proj, roi_flat, toc_bytes )
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
% toc_bytes : struct with data transfers during parpool jobs
%
% RETURNS
% proj : flat-field corrected projections
% corr : struct with correlation information
% toc_bytes : struct with data transfers during parpool jobs
%
% Written by Julian Moosmann.
%
% [proj, corr] = proj_flat_correlation( proj, flat, image_correlation, par, roi_proj, roi_flat )

%% To do
% Combine structs found at different locations and check what happens when
% different

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    roi_proj = [];
end
if nargin < 6
    roi_flat = [];
end

%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters in struct fields
force_calc = image_correlation.force_calc;
method = image_correlation.method;
raw_bin = par.raw_bin;
flat_corr_area1 = image_correlation.area_width;
flat_corr_area2 = image_correlation.area_height;
corr_num_flats = image_correlation.num_flats;
[im_shape_cropbin1, im_shape_binned2, num_proj_used] = size( proj );
[im_shape_binned1,~,num_ref_used] = size( flat );
x0 = par.offset_shift_x0;
%x1 = par.offset_shift_x1;
outpath = write.parpath;

%% CORRELATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
    
    case {'', 'median'}
        fprintf( 'No correlation.' )
        corr = [];
        
    otherwise
        
        % Create corr struct from current parameters
        corr_curr.datetime = datetime;
        corr_curr.par.method = method;
        corr_curr.par.raw_roi = par.raw_roi;
        corr_curr.par.raw_bin = raw_bin;
        corr_curr.par.ref_range = par.ref_range;
        corr_curr.par.proj_range = par.proj_range;
        corr_curr.par.size.proj = size( proj );
        corr_curr.par.size.flat = size( flat );
        corr_curr.par.flat_corr_area1 = flat_corr_area1;
        corr_curr.par.flat_corr_area2 = flat_corr_area2;
        corr_curr.par.image_correlation_filter = image_correlation.filter;
        if image_correlation.filter
            corr_curr.par.image_correlation_filter_type = image_correlation.filter_type;
            corr_curr.par.image_correlation_filter_parameter = image_correlation.filter_parameter;
        end
        
        %% Search for previously calculated correlation matrix
        % Loop over processed and scratch_cc
        match = 0;
        filename = sprintf( '%scorr.mat', outpath);
        if regexp( 'processed', filename )
            filename_proc = filename;
            filename_scra = regexprep( filename, 'processed', 'scratch_cc' );
        else
            filename_proc = regexprep( filename, 'scratch_cc', 'processed' );
            filename_scra = filename;
        end
        
        % Check for file in both directories
        for filename = { filename_proc, filename_scra }
            
            % Load correlation matrix if exists, otherwise break iteration
            if exist( filename{1} , 'file' )
                load( filename{1},  'corr' );
            else
                continue
            end
            
            % Check reco paramaters
            if exist( 'corr', 'var' )
                
                % Loop over structs
                for nn = 1:numel( corr )
                    
                    % Compare current and found parameters
                    match = isequal( corr(nn).par, corr_curr.par );
                    %fprintf( '\n file: %s, struc_ind: %u, match: %u', filename{1}, nn, match )
                    
                    % exit loop over found structs if match
                    if match
                        match_ind = nn;
                        %fprintf( '\n Match: exit struct loop ' )
                        break
                    end
                end % for nn = 1:numel( corr )
            end % if exist( 'corr', 'var' )
            
            % Exit file loop if match
            if match
                %fprintf( '\n Match: Exit file loop' )
                break
            end
        end
        
        % Correlation matrix found -> run correlation not necessary
        run_corr = ~match;
        
        % Force run correlation
        if force_calc
            run_corr = 1;
        end
        
        % Print info
        if match && ~force_calc
            fprintf( '\nLoading correlation computed on %s!', corr(match_ind).datetime )
            corr_mat = corr(match_ind).mat;
        end
        
        %% Correlate
        if run_corr
            
            % Correlation ROI
            flat_corr_area1 = IndexParameterToRange( flat_corr_area1, im_shape_binned1 );
            flat_corr_area2 = IndexParameterToRange( flat_corr_area2, im_shape_binned2 );
            
            % ROI
            if isempty( roi_proj ) && isempty( roi_flat )
                roi_proj = proj(flat_corr_area1,flat_corr_area2,:);
                roi_flat = flat(flat_corr_area1,flat_corr_area2,:);
            end
            
            % Correlate flat fields: Compute correlation for each pair projection/flat-field
            fprintf( '\nCorrelate projections & flat-fields. Method: %s.', method )
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
            
            % Add calculated matrix to current struct
            corr_curr.mat = corr_mat;
            
            % Append current correlation struct to found struct
            if ~exist( 'corr', 'var' )
                corr(1) = corr_curr;
            else
                if match
                    corr(match_ind) = corr_curr;
                else
                    corr(numel( corr ) + 1) = corr_curr;
                end
            end
            
            % Save correlation matrix
            CheckAndMakePath( outpath )
            filename = sprintf( '%scorr.mat', outpath);
            save( filename, 'corr' )
            
            fprintf( ' Done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
            fprintf( '\n correlation matrix path :\n  %s', filename )
            
        end %  if ~force_calc
        
end % switch method

%% Plot correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.visual_output && exist( 'corr_mat', 'var' )
    figure( 'Name', 'Correlation of projections and flat-fields');
    mid = round( num_proj_used / 2 );
    f = @(mat,pp) mat(pp,:)';
    Y = [ f( corr_mat, 1), f( corr_mat, mid), f( corr_mat, num_proj_used) ];
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
    case 'none'
        fprintf( '\nNO FLAT FIELD CORRECTION')
        
    case {'', 'median'}
        
        % Flat field correction without correlation
        fprintf( '\nFlat-field correction w/o correlation.')
        
        flat_median = median( flat, 3);
        parfor nn = 1:num_proj_used
        %    im = proj(:, :, nn);
            %flat_median_shifted = circshift( flat_median, round( -x0(nn) / raw_bin ), 1 );
            %proj(:, :, nn) = im ./ flat_median_shifted;
            
            % Binned shift (shift, not first pixel)
            shift = ( x0(nn) - 1 ) / raw_bin;
            shift_int = floor( shift );
            shift_sub = shift - shift_int;
            
            % shift flat
            if mod( shift_sub, 1 ) ~= 0
                % crop flat at integer shift, then shift subpixel
                xx = shift_int + (1:im_shape_cropbin1+1);
                flat_median_shifted = imtranslate( flat_median(xx,:), [0 -shift_sub], 'linear' );
            else
                xx = shift_int + (1:im_shape_cropbin1);
                flat_median_shifted = flat_median(xx,:);
            end
            
            % flat field correction
            p = proj(:, :, nn);
            p = p ./ flat_median_shifted(1:im_shape_cropbin1,:) ;
            % Reassign
            proj(:,:,nn) = p;
            
        end
        
    otherwise
        fprintf( '\nFlat-field correction using %u best match(es).', corr_num_flats)
        t = toc;
        
        [~, corr_mat_pos] = sort( normat( corr_mat ), 2);
        flat_ind = corr_mat_pos(:,1:corr_num_flats);
        
        %        startS = ticBytes( gcp );
%         [~, mem_avail_cpu, ~] = free_memory;
%         mem_proj = Bytes( proj );
%         mem_flat = Bytes( flat );
%         parpool_max = min( [floor( mem_avail_cpu / ( mem_proj + mem_flat ) ), par.poolsize] );
%         parfor ( nn = 1:num_proj_used, parpool_max )
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
        %toc_bytes.correction = tocBytes( gcp, startS );
        
        fprintf( ' Done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
