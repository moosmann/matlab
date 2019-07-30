function [proj, corr, roi] = proj_flat_correlation( proj, flat, image_correlation  )
% Correlate all projection with all flat-field to find the best matching
% pairs for the flat-field correction using method of choice.
%
% proj : stack of projections
% flat : stack of flat field
% image_correlation : struct, required fields see below
%
% Written by Julian Moosmann.
%
% proj_flat_correlation( proj, flat, image_correlation  )

%% TODO: Clean up and simplify

% Parameters in struct fields
method = image_correlation.method;
flat_corr_area1 = image_correlation.area_width;
flat_corr_area2 = image_correlation.area_height;
corr_shift_max_pixelshift = image_correlation.shift.max_pixelshift;
corr_num_flats = image_correlation.num_flats;
decimal_round_precision = image_correlation.decimal_round_precision;
force_calc = image_correlation.force_calc;
im_shape_binned = image_correlation.im_shape_binned;
flatcor_path = image_correlation.flatcor_path;
verbose = image_correlation.verbose;
visual_output = image_correlation.visual_output;
poolsize = image_correlation.poolsize;

num_proj_used = size( proj, 3);
num_ref_used = size( flat, 3);

%% Projection / flat-field correlcation
switch method
    
    case {'none', ''}
        % Flat field correction without correlation
        flat_median = median( flat, 3);
        proj = bsxfun( @times, proj, 1./flat_median);
        corr = [];
        roi_proj =[];
        roi_flat = [];
        
    otherwise
        % Load previously calculated correlation matrix & check if valid
        if ~force_calc
            filename_scra = sprintf( '%scorr.mat', flatcor_path);
            filename_proc = regexprep( filename_scra, 'scratch_cc', 'processed' );
            if exist( filename_proc, 'file' )
                load( filename_proc,  'corr' );
            else
                if exist( filename_scra, 'file' )
                    load( filename_scra,  'corr' );
                end
            end
            if exist( 'corr', 'var' )
                force_calc = ~(isequal( corr.method, method) ...
                    && isequal( corr.im_shape_binned, im_shape_binned ) ...
                    && isequal( corr.size.proj, size( proj ) ) ...
                    && isequal( corr.size.flat, size( flat ) ) ) ...
                    && isequal( corr.raw_roi, image_correlation.raw_roi ) ...
                    && isequal( corr.raw_bin, image_correlation.raw_bin ) ...
                    && isequal( corr.bin_before_filtering, image_correlation.bin_before_filtering ) ...
                    && isequal( corr.proj_range, image_correlation.proj_range ) ...
                    && isequal( corr.ref_range, image_correlation.ref_range );
            else
                force_calc = 1;
            end
        end
        
        % Correlation ROI
        flat_corr_area1 = IndexParameterToRange(flat_corr_area1, im_shape_binned(1));
        flat_corr_area2 = IndexParameterToRange(flat_corr_area2, im_shape_binned(2));
        
        roi_flat = flat(flat_corr_area1, flat_corr_area2, :);
        roi_proj = proj(flat_corr_area1, flat_corr_area2, :);
        
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
                f = roi_flat(:,:,ff);
                
                parfor pp = 1:num_proj_used
                    
                    % projection
                    p = roi_proj(:,:,pp);
                    
                    if sum(strcmpi( method, {'shift','all'}))
                        % shift via image cross correlation
                        out = ImageCorrelation( p, f, 0, 0, 0, 0, 1);
                        c_shift_1(pp,ff) = round( out.shift1, decimal_round_precision );
                        c_shift_2(pp,ff) = round( out.shift2, decimal_round_precision) ; % relevant shift
                        
                    elseif sum(strcmpi( method, {'diff','all'}))
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
                        
                    elseif sum(strcmpi( method, {'cov', 'corr', 'ssim','all'}))
                        % input
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
                        c_ssim_ml(pp,ff) = - ssim( roi_proj(:,:,pp), f ); %'DynamicRange', 'uint16'
                    end
                end
                
                switch method
                    case 'shift'
                        corr.mat = c_shift_2;
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
                    case 'ssim-ml'
                        corr.mat = c_ssim_ml;
                    case 'cov'
                        corr.mat = c_cov;
                    case 'corr'
                        corr.mat = c_corr;
                    case 'all'
                        corr.shift1 = c_shift_1;
                        corr.shift2 = c_shift_2;
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
                corr.im_shape_binned = im_shape_binned;
                corr.size.proj = size( proj );
                corr.size.flat = size( flat );
                corr.datetime = datetime;
                corr.raw_roi = image_correlation.raw_roi;
                corr.raw_bin = image_correlation.raw_bin;
                corr.bin_before_filtering = image_correlation.bin_before_filtering;
                corr.proj_range = image_correlation.proj_range;
                corr.ref_range = image_correlation.ref_range;
                CheckAndMakePath( flatcor_path )
                save( sprintf( '%scorr.mat', flatcor_path), 'corr' )
                
            end
            PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        end
        
        %% Visual output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if visual_output
            figure('Name', 'Correlation of projections and flat-fields');
            
            switch method
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
        
        %% Flat- and dark-field correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            case 'shift'
                % best match
                if corr_shift_max_pixelshift == 0
                    [~, pos] = min( abs( corr.mat ), [], 2 );
                    
                    if pflag
                        parfor nn = 1:num_proj_used
                            f = flat(:, :, pos(nn));
                            proj(:, :, nn) = proj(:, :, nn) ./ f;
                        end
                    else
                        for nn = 1:num_proj_used
                            f = flat(:, :, pos(nn));
                            proj(:, :, nn) = proj(:, :, nn) ./ f;
                        end
                    end
                    
                    % use all flats which are shifted less pixels than corr_shift_max_pixelshift
                elseif corr_shift_max_pixelshift > 0
                    nflats = zeros(1, num_proj_used);
                    if pflag
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
                            f = squeeze( mean( flat(:, :, flat_ind), 3) );
                            proj(:, :, nn) = proj(:, :, nn) ./ f;
                        end
                    else
                        for nn = 1:num_proj_used
                            vec = 1:num_ref_used;
                            flat_ind = vec( abs( c_shift_2(nn, :) ) < corr_shift_max_pixelshift );
                            if numel( flat_ind ) > corr_num_flats
                                flat_ind( corr_num_flats + 1:end ) = [];
                            end
                            if isempty( flat_ind )
                                [~, flat_ind] = min( abs( c_shift_2(nn, :) ) );
                            end
                            nflats(nn) = numel(flat_ind);
                            f = squeeze( mean( flat(:, :, flat_ind), 3) );
                            proj(:, :, nn) = proj(:, :, nn) ./ f;
                        end
                    end
                    
                    PrintVerbose(verbose, '\n number of flats used per projection: [mean, min, max] = [%g, %g, %g]', mean( nflats ), min( nflats ), max( nflats) )
                else
                    error('Value of maximum shift (%g) is not >= 0', corr_shift_max_pixelshift)
                end
                
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
                    parfor nn = 1:num_proj_used
                        f = squeeze( mean( flat(:, :, flat_ind(nn,:)), 3) );
                        proj(:, :, nn) = proj(:, :, nn) ./ f;
                    end
                else
                    for nn = 1:num_proj_used
                        f = squeeze( mean( flat(:, :, flat_ind(nn,:)), 3) );
                        proj(:, :, nn) = proj(:, :, nn) ./ f;
                    end
                end
                
                %         % Switch back to old poolsize
                %         if num_worker < poolsize
                %             fprintf( '\n Switch back to old pool size %u -> %u.', num_worker, poolsize )
                %             OpenParpool( poolsize, 0, [beamtime_path filesep 'scratch_cc']);
                %         end
                
                PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        end
end
roi.proj = roi_proj;
roi.flat = roi_flat;
