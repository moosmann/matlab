function [proj, corr] = proj_flat_correlation( proj, flat, image_correlation, par, roi_proj  )
% Correlate all projection with all flat-field to find the best matching
% pairs for the flat-field correction using method of choice.
%
% proj : stack of projections
% flat : stack of flat field
% image_correlation : struct, required fields see below
% par : struct, required fields see below
% roi_proj : ROI for correlation when offset shift is corrected
%
% Written by Julian Moosmann.
%
% proj_flat_correlation( proj, flat, image_correlation, par, roi_proj )

if nargin < 5
    roi_proj = [];
end

%% TODO: Clean up and simplify

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
visual_output = par.visual_output;
poolsize = par.poolsize;
x0 = par.offset_shift_x0;
%x1 = par.offset_shift_x1;
raw_bin = par.raw_bin;

im_shape_cropbin1 = size( proj, 1 );
num_proj_used = size( proj, 3);
num_ref_used = size( flat, 3);

%% CORRELATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
    
    case {'none', '', 'median'}
        % Flat field correction without correlation
        flat_median = median( flat, 3);
        parfor nn = 1:num_proj_used
            im = proj(:, :, nn);
            flat_median_shifted = circshift( flat_median, round( -x0(nn) / raw_bin ), 1 );
            proj(:, :, nn) = im ./ flat_median_shifted;
        end
        corr = [];
        
    otherwise
        % Load previously calculated correlation matrix & check if valid
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
        
        % Correlation ROI
        flat_corr_area1 = IndexParameterToRange(flat_corr_area1, im_shape_binned1 );
        flat_corr_area2 = IndexParameterToRange(flat_corr_area2, im_shape_binned2 );
        
        roi_flat = flat(flat_corr_area1, flat_corr_area2, :);
        if isempty( roi_proj )
            roi_proj = proj(flat_corr_area1, flat_corr_area2, :);
        end
        
        if ~force_calc
            dt = '??';
            if isfield( corr, 'datetime' )
                dt = corr.datetime;
            end
            PrintVerbose( verbose, '\nLoading correlation computed on %s!', dt )
        else
            
            % Correlate flat fields
            PrintVerbose(verbose, '\nCorrelate projections and flat-fields. Method: %s.', method)
            t = toc;
            
            % Preallocation
            foo = zeros( num_proj_used, num_ref_used);
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
            
            % Compute correlation for each pair projection/flat-field
            
            if strcmpi( method, 'ssim-g' )
                radius = 1.5;
                filtRadius = ceil(radius*3); % 3 Standard deviations include >99% of the area.
                filtSize = 2*filtRadius + 1;
                gaussFilterFcn = @(X)imgaussfilt(X, radius, 'FilterSize', filtSize, 'Padding','replicate');
            end
            
            for ff = 1:num_ref_used
                
                % flat field
                f = roi_flat(:,:,ff);
                
                if  sum(strcmpi( method, {'cov', 'corr', 'ssim','ssim-g', 'all'}))
                    f_mean = mean2( f );
                    f_std = std2( f );
                end
                if strcmpi( method, 'ssim-g' )
                    f = gaussFilterFcn( f );
                end
                
                %parfor pp = 1:num_proj_used
                for pp = 1:num_proj_used
                    
                    % projection
                    p = roi_proj(:,:,pp);
                    
                    if sum(strcmpi( method, {'diff','all'}))
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
                        
                    elseif sum(strcmpi( method, {'std','all'}))
                        % std : standard deviation of ratio of p and f
                        c_std(pp,ff) = std2( p ./ f );
                        
                    elseif sum(strcmpi( method, {'entropy','all'}))
                        % entropy : entropy of ratio of p and f
                        c_ent(pp,ff) = entropy( double( p ./ f ) );
                        
                    elseif sum(strcmpi( method, {'cross-entropy12', 'cross-entropy21', 'cross-entropyx','all'}))
                        % cross entropy of proj and flat
                        p1 = imhist( normat( p(:) ) ) ./ numel( p );
                        p2 = imhist( normat( f(:) ), numel( p1 ) ) ./ numel( f );
                        m = boolean( (p1 == 0) + (p2 == 0) );
                        p1(m) = [];
                        p2(m) = [];
                        c_cross_entropy12(pp,ff) = - sum( p1 .* log( p2 ) );
                        c_cross_entropy21(pp,ff) = - sum( p2 .* log( p1 ) );
                        c_cross_entropyx(pp,ff) = sum( p2 .* log2( p1 ) - p1 .* log2( p2 ) );
                        
                    elseif sum(strcmpi( method, {'cov', 'corr', 'ssim','ssim-g','all'}))
                        % input
                        if strcmpi( method, 'ssim-g' )
                            p = gaussFilterFcn( p );
                        end
                        p_mean = mean2( p );
                        p_std = std2( p );
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
                        % L = round(max(max(roi_flat(:)),max(roi_flat(:))) - min(min(roi_flat(:)),min(roi_flat(:))));
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
                        
                    elseif sum(strcmpi( method, {'ssim-ml','all'}))
                        % Matlab's structural similarity index (SSIM)
                        c_ssim_ml(pp,ff) = - ssim( p, f ); %'DynamicRange', 'uint16'
                    end
                end
                
                switch method
                    case 'diff'
                        corr.mat = c_diff1_l2;
                    case 'std'
                        corr.mat = c_std;
                    case 'entropy'
                        corr.mat = c_ent;
                    case 'cross-entropy12'
                        corr.mat = c_cross_entropy12;
                    case 'cross-entropy21'
                        corr.mat = c_cross_entropy21;
                    case 'cross-entropyx'
                        corr.mat = c_cross_entropyx;
                    case 'ssim'
                        corr.mat = c_ssim;
                    case 'ssim-g'
                        corr.mat = c_ssim;
                    case 'ssim-ml'
                        corr.mat = c_ssim_ml;
                    case 'cov'
                        corr.mat = c_cov;
                    case 'corr'
                        corr.mat = c_corr;
                    case 'all'
                        corr.diff1_l1 = c_diff1_l1;
                        corr.diff1_l2 = c_diff1_l2;
                        corr.diff2_l1 = c_diff2_l1;
                        corr.diff2_l2 = c_diff2_l2;
                        corr.std = c_std;
                        corr.entropy = c_ent;
                        corr.cross_entropy12 = c_cross_entropy12;
                        corr.cross_entropy21 = c_cross_entropy21;
                        corr.cross_entropyx = c_cross_entropyx;
                        corr.ssim = c_ssim;
                        corr.ssim_ml = c_ssim_ml;
                        corr.cov = c_cov;
                        corr.corr = c_corr;
                end
                
                if strcmp( method, 'all' )
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
                % Save correlation matrix
                %corr.mat = corr_mat;
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
                
            end
            PrintVerbose(verbose, ' Done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        end
        
        %% Visual output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if visual_output
            figure('Name', 'Correlation of projections and flat-fields');
            
            switch method
                case {'cov', 'corr', 'ssim', 'ssim-g'}
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
                    %% TODO
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
                    Y = [ f( corr.mat, 1), f( corr.mat, mid), f( corr.mat, num_proj_used) ];
                    plot( Y, '.' )
                    legend( 'first proj', 'mid proj', 'last proj' )
                    axis tight
                    title(sprintf('correlation method: %s', method))
                    xlabel( 'flat field index' )
                    ylabel( 'measure' )
            end
            drawnow
        end
        
        %% CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        PrintVerbose(verbose, '\nFlat- and dark-field correction.')
        t = toc;
        
        % Memory requirement and parpool processing
        [mem_free, ~, mem_total] = free_memory;
        mem_req_tot = Bytes( proj ) + Bytes(flat) * (poolsize + 1);
        mem_req_par = Bytes(flat) * poolsize;
        if mem_req_tot < 0.75 * mem_total && mem_req_par < mem_free
            pflag = 1;
            fprintf( ' Use parpool.' )
        else
            pflag = 0;
            fprintf( ' No parpool.' )
        end
        
        % Reduce number of workers because of memory overhead
        %         num_worker =  max( floor( 0.5 * min( mem_free/1024^3, mem_total / 1024^3 - GB(proj) - GB(flat) ) / ceil( GB(flat) ) ) - 1, 1 );
        %         if num_worker < poolsize
        %             fprintf( '\n Change pool size %u -> %u for flat field correction. (Each worker requires a full copy of flat fields).', poolsize, num_worker)
        %             OpenParpool( num_worker, 0, [beamtime_path filesep 'scratch_cc'], 1);
        %         end
        
        switch method
            case 'all'
                % Save correlation matrix
                CheckAndMakePath( flatcor_path )
                save( sprintf( '%s/corr_all.mat', flatcor_path), 'corr' )
                fprintf( 'All available measures calulated. NO FLAT FIELD CORRECTION DONE.')
                
            otherwise
                [~, corr_mat_pos] = sort( normat( corr.mat ), 2);
                % Flat field correction
                flat_ind = corr_mat_pos(:,1:corr_num_flats);
                if pflag
                    %parfor nn = 1:num_proj_used
                    for nn = 1:num_proj_used
                        % Mean over flat fields
                        flat_mean = squeeze( mean( flat(:, :, flat_ind(nn,:)), 3) );
                        % Binned shift (shift not first pixel)
                        shift = (x0(nn) - 1 ) / raw_bin;
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
                else
                    for nn = 1:num_proj_used
                        % Mean over flat fields
                        flat_mean = squeeze( mean( flat(:, :, flat_ind(nn,:)), 3) );
                        % Binned shift
                        shift = (x0(nn) - 1 ) / raw_bin;
                        shift_int = floor( shift );
                        shift_sub = shift - shift_int;
                        % the subpixel offest is tranlate back!
                        shift_sub = shift_sub - 1;
                        
                        % shift flat
                        if mod( shift_sub, 1 ) ~= 0
                            % crop flat at integer shift, then shift subpixel
                            xx = shift_int + (1:im_shape_cropbin1+1);
                            flat_mean_shifted = imtranslate( flat_mean(xx,:), [0 shift_sub], 'linear' );
                        else
                            xx = shift_int + (1:im_shape_cropbin1);
                            flat_mean_shifted = flat_mean(xx,:);
                        end
                        
                        % flat field correction
                        p = proj(:, :, nn);
                        p = p ./ flat_mean_shifted(1:im_shape_cropbin1,:) ;
                        proj(:,:,nn) = p;
                    end
                end
                
                %         % Switch back to old poolsize
                %         if num_worker < poolsize
                %             fprintf( '\n Switch back to old pool size %u -> %u.', num_worker, poolsize )
                %             OpenParpool( poolsize, 0, [beamtime_path filesep 'scratch_cc']);
                %         end
                
                PrintVerbose(verbose, ' Done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        end
end
