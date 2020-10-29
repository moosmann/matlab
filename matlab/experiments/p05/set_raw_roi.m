function par = set_raw_roi( par, im_raw,  scan_path, fig_path, ref_full_path, dark_names )
% Set raw ROI.
%
% ARGUMENT
% raw_roi : default: [], defines the vertical (and horizontal) range of the
%   image to be read in. Supports absolute, relative and negative indexing.
%   Negative indexing only for the 2nd/end index in each dimension.
%   Horizontal ROI is only supporte for a TIFF images if at all.
%   if [] : use full image
%   if = -1 : Auto ROI for DCM only. Values < -1 increases the threshold.
%   if raw_roi >= 1 : interpreted as absolute values.
%   if raw_roi > 0 and < 1 : interpreted as relative values w.r.t. to the
%      image height
%   if raw_roi(2) or raw_roi(4) < 0 : values are counted from the end of
%      the index range.
%   if raw_roi(2) or raw_roi(4) > -1 and < 0 : negative relative indexing.
%   if raw_roi(2) or raw_roi(4) <= -1 : absolute negative indexing couning
%      from the end 
%   
% Written by Julian Moosmann.
%
% raw_roi = set_raw_roi( raw_roi, par, im_shape_raw, im_raw, tif_info, dtype, im_trafo, scan_path, fig_path, ref_names, dark_names )
% raw_roi = set_raw_roi( raw_roi, par, im_shape_raw, im_raw, tif_info, dtype, im_trafo, scan_path, fig_path, ref_names, dark_names )

%% Default arguments
assign_from_struct( par, 'raw_roi',  []);
assign_from_struct( par, 'im_shape_raw',  []);
assign_from_struct( par, 'tif_info',  []);
assign_from_struct( par, 'dtype',  []);
assign_from_struct( par, 'im_trafo',  []);


% old: par.raw_roi = set_raw_roi( par.raw_roi, par, par.im_shape_raw, im_raw, par.tif_info, par.dtype, par.im_trafo, scan_path, fig_path, ref_names, dark_names );
par.raw_roi_old = raw_roi;
if ~isempty( raw_roi ) % else AUTO ROI
    if numel( raw_roi ) > 1
        % relative indexing
        for nn = numel( raw_roi ):-1:1
            x = raw_roi(nn);
            if x >= -1 && x <= 1
                switch nn
                    case 1
                        if x < 1
                            x = round( ( im_shape_raw(2) - 1) * x + 1, -1 ) + 1;
                        end
                    case 2
                        if x < 0
                            x = 1 + x;
                        end
                        x = round( ( im_shape_raw(2) - 1) * x + 1, -1 );
                    case 3
                        if x < 1
                            x = round( ( im_shape_raw(1) - 1) * x + 1, -1 ) + 1;
                        end
                    case 4
                        if x < 0
                            x = 1 + x;
                        end
                        x = round( ( im_shape_raw(1) - 1) * x + 1, -1 );
                end
            end
            raw_roi(nn) = x;
        end
        
        if raw_roi(2) <= 0 % negative index counting from the end
            raw_roi(2) = im_shape_raw(2) + raw_roi(2);
        end
        
        if numel( raw_roi ) == 4
            if raw_roi(4) <= 0 % negative index counting from the end
                raw_roi(4) = im_shape_raw(1) + raw_roi(4);
            end
        end
    else
        if raw_roi(1) < 1 % indicate AUTO ROI
            % Make roi_fac dependent on dark field
            if raw_roi(1) == 0
                roi_fac = 4;
            else
                roi_fac = abs( raw_roi );
            end
            % Read first non-zero flat
            mm = 1;
            while mean2( im_raw ) == 0
                mm = mm + 1;
                %filename = sprintf('%s%s', scan_path, ref_names{mm});
                filename = ref_full_path{mm};
                %im_raw = read_image( filename, '', [], tif_info, im_shape_raw, dtype, im_trafo );
                im_raw = read_image( filename, par, 1 );
            end
            % Read last non-zero flat
            mm = numel( ref_full_path );
            %filename = sprintf('%s%s', scan_path, ref_names{mm});
            filename = ref_full_path{mm};
            %im_raw2 = read_image( filename, '', [], tif_info, im_shape_raw, dtype, im_trafo );
            im_raw2 = read_image( filename, par, 1 );
            while mean2( im_raw2 ) == 0
                mm = mm - 1;
                %filename = sprintf('%s%s', scan_path, ref_names{mm});
                filename = ref_full_path{mm};
                %im_raw2 = read_image( filename, '', [], tif_info, im_shape_raw, dtype, im_trafo );
                im_raw2 = read_image( filename, par, 1);
            end
            % Read non-zero dark
            mm = 1;
            im_raw_dark = 0;
            while mean2( im_raw_dark ) == 0
                mm = mm + 1;
                filename = sprintf('%s%s', scan_path, dark_names{mm});
                %im_raw_dark = read_image( filename, '', [], tif_info, im_shape_raw, dtype, im_trafo );
                im_raw_dark = read_image( filename, par, 1 );
            end
            % Threshold
            im_raw_line = median( im_raw / 2 + im_raw2 / 2, 1);
            roi_thresh = double( min( im_raw_line ) + roi_fac * median( im_raw_dark(:) ) );
            % Crop indices
            pl = 1;
            while im_raw_line(pl) < roi_thresh && pl < im_shape_raw(2) / 2 - 100
                pl = pl + 1;
            end
            pl = 10 * ceil( pl / 10 );
            pr = im_shape_raw(2);
            while im_raw_line(pr) < roi_thresh && pr > im_shape_raw(2) / 2 + 100
                pr = pr - 1;
            end
            pr = 10 * floor( pr / 10 );
            raw_roi = double( [pl pr] );
            fprintf( '\n vertical auto roi : [%u %u], roi threshold factor : %g', raw_roi, roi_fac )
            
            % Plot auto ROI
            if par.visual_output
                f = figure( 'Name', 'auto ROI: raw image and cropping region' );
                
                subplot(1,2,1)
                im = (single(im_raw) + single(im_raw2) ) / 2;
                pll = max( 1, pl - 20 );
                prr = min( im_shape_raw(2), pr + 20);
                im(:,[pll:pl, pr:prr] ) = max( im(:) );
                imsc(im);
                title(sprintf('raw flat field: first + last'))
                axis equal tight
                xticks('auto'),yticks('auto')
                camroll(90)
                
                subplot(1,2,2)
                plot( [im_raw_line' repmat( roi_thresh, [numel(im_raw_line) 1]) ] )
                camroll(90)
                title(sprintf('horizontal projection and cut level\nthreshold factor: %g', roi_fac ))
                axis tight
                ax = gca;
                ax.YAxisLocation = 'right';
                ax.XDir = 'reverse';
                ax.YDir = 'normal';
                text( raw_roi(1), roi_thresh + 80, sprintf('raw roi(1)=%u', raw_roi(1) ) )
                text( raw_roi(2), roi_thresh + 80, sprintf('raw roi(2)=%u', raw_roi(2) ) )
                
                drawnow
                saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
            end % if par.visual_output Plot auto ROI
        end
    end
end

par.raw_roi = raw_roi;
