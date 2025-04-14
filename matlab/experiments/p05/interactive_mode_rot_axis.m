function [tomo, angles, tint, par] = interactive_mode_rot_axis( par, logpar, phase_retrieval, tomo, write, interactive_mode, proj, angles)
% Interactive mode to determine the position and tilt (optionally) of the
% rotation axis. When position and tilt are searched for, two intertwined
% loops are used alternating stiltthe determination of the position and the tilt
% of the rotation axis.
%
% Written by Julian Moosmann
%
% [tomo, tint] = interactive_mode_rot_axis( par, logpar, phase_retrieval, tomo, write, interactive_mode, proj, angles)

tomo.angle_scaling = 1;
%tomo.angles = angles;
imsc1 = @(im) imsc( rot90( im ) );
tint = 0;
angle_scaling = [];
lamino = interactive_mode.lamino;
tomo.lamino = lamino;
show_stack_imagej = interactive_mode.show_stack_imagej;
show_stack_imagej_use_virtual = interactive_mode.show_stack_imagej_use_virtual;
window_state = par.window_state;
if isempty( tomo.rot_axis_offset )
    tomo.rot_axis_offset = 0;
end

if tomo.run || tomo.run_interactive_mode
    
    fprintf( '\nTomography:')
    
    % Full rotation angle
    if isempty( tomo.rot_angle_full_range )
        if ~isempty( logpar )
            if isfield( logpar, 'rotation')
                % From log file
                tomo.rot_angle_full_range = logpar.rotation / 180 * pi;
            elseif exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
                % from beam current log
                tomo.rot_angle_full_range = (cur.proj(end).angle - cur.proj(1).angle) * pi /180; % KIT: , EHD: ok
            end
        end
    end
    if isempty( tomo.rot_angle_full_range )
        tomo.rot_angle_full_range = max( angles(:) ) - min( angles(:) );
    end
    if isscalar( tomo.rot_angle_full_range)
        ar = tomo.rot_angle_full_range;
    else
        ar = max(tomo.rot_angle_full_range) - min(tomo.rot_angle_full_range);
    end
    fprintf('\n rotation angle range: %f * pi = %f deg', ar / pi, ar / pi * 180)
    
    if numel( angles ) ~= size( proj, 3 )
        error('Number of elements in array of angles (%g) unequal number of projections read (%g)', numel( angles ), size( proj, 3))
    end
    
    % retrieve index at angles 0 and pi

    [~, ind1] = min( abs( angles ));
    [~, ind2] = min( abs( angles - pi));
        
    [im_shape_cropbin1, im_shape_binned2, ~] = size( proj );
    
    tomo.rot_axis_position = im_shape_cropbin1 / 2 + tomo.rot_axis_offset;
    
    % Tilt of rotation axis
    if interactive_mode.rot_axis_tilt
        
        % ROI for correlation of projections at angles 0 & pi
        if isempty( tomo.rot_axis_corr_area1 )
            r = tomo.rot_axis_position / im_shape_cropbin1;
            if r < 1 / 4
                tomo.rot_axis_corr_area1 = [0 0.25];
            elseif r > 3 / 4
                tomo.rot_axis_corr_area1 = [0.75 1];
            else
                tomo.rot_axis_corr_area1 = [0.25 0.75];
            end
        end
        tomo.rot_axis_corr_area1 = IndexParameterToRange(tomo.rot_axis_corr_area1,im_shape_cropbin1);
        tomo.rot_axis_corr_area2 = IndexParameterToRange(tomo.rot_axis_corr_area2,im_shape_binned2);        
        if max(angles) > 1.9 * pi
            im1c = proj(:,tomo.rot_axis_corr_area2,ind1);
            im2c = proj(:,tomo.rot_axis_corr_area2,ind2);
        else
            im1c = RotAxisSymmetricCropping(proj(:,tomo.rot_axis_corr_area2,ind1),tomo.rot_axis_position,1);
            im2c = flipud(RotAxisSymmetricCropping(proj(:,tomo.rot_axis_corr_area2,ind2),tomo.rot_axis_position,1));
        end
        [optimizer, metric] = imregconfig('monomodal');
        tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);
        rot_axis_tilt_calc = asin( tform_calc.T(1,2) ) / 2;
    else
        rot_axis_tilt_calc = [];
    end
    
    %% INTERACTIVE MODE: rotation axis position / tilt %%%%%%%%%%%%%%%%%%%%%
    
    tint = 0;
    if interactive_mode.rot_axis_pos
        tint = toc;
        cprintf('RED','\nSTART OF INTERACTIVE MODE')
        
        % parameter strcut for interactive functions
        itomo = tomo;
        if isempty( tomo.take_neg_log )
            if phase_retrieval.apply && phase_retrieval.apply_before
                itomo.take_neg_log = 0;
            else
                itomo.take_neg_log = 1;
            end
        end
        if phase_retrieval.apply
            itomo.inumber_of_stds = 9;
        else
            itomo.inumber_of_stds = 4;
        end
        if interactive_mode.slice_number > 1 && interactive_mode.slice_number <= size( proj, 2 )
            slice = interactive_mode.slice_number;
        elseif interactive_mode.slice_number <= 1 && interactive_mode.slice_number >= 0
            slice = floor((size( proj, 2 ) - 1) * interactive_mode.slice_number + 1 );
        end
        
        itomo.offset = tomo.rot_axis_offset;
        if ~lamino
            itomo.tilt = tomo.rot_axis_tilt_camera;
            itomo.fixed_tilt = tomo.rot_axis_tilt_lamino;
        else
            itomo.tilt = tomo.rot_axis_tilt_lamino;
            itomo.fixed_tilt = tomo.rot_axis_tilt_camera;
        end
        
        itomo.slice = slice;
        % Initialize paramters to be determined in interactive loops
        offset = [];
        tilt = [];
        angle_scaling = [];
        
        % Loop over offsets
        while ~isscalar( offset )
            interactive_mode.rot_axis_search_range = offset;
            cprintf( 'RED', '\nEntering offset loop' )
            % Print parameters
            fprintf('\n scan: %s', write.scan_name )
            fprintf('\n number of pixels: %u', im_shape_cropbin1)
            fprintf('\n image center: %.1f', im_shape_cropbin1 / 2)
            fprintf('\n current slice : %u', slice)
            fprintf('\n current rotation axis position : %.2f', itomo.rot_axis_position)
            fprintf('\n current rotation axis TILT : %g rad = %g deg', itomo.tilt, itomo.tilt * 180 / pi)
            fprintf('\n current rotation axis OFFSET : %.2f', itomo.rot_axis_offset )
            fprintf('\n default offset range : current OFFSET + [')
            fprintf(' %.2g', interactive_mode.rot_axis_pos_default_search_range )
            fprintf(']' )
            
            % Query parameters
            inp = '';
            while ischar( inp )
                switch inp
                    case 's'
                        slice = input( sprintf( '\nENTER ABSOLUTE [1,%u] OR RELATIVE [0,1] SLICE NUMBER : ', size( proj, 2 )) );
                        if slice <= 1 && slice >= 0
                            slice = round((size( proj, 2 ) - 1) * slice + 1 );
                        end
                        fprintf( ' new slice : %u (before: %u)', slice, itomo.slice );
                        itomo.slice = slice;
                        tomo.slice = slice;
                    case 'd'
                        cprintf( 'RED', 'ENTERING DEBUG MODE. Continue with F5.' )
                        keyboard
                    case 'c'
                        txt = sprintf( '\nENTER DISTANCE TO OUTER OFFSET:');
                        tmp = input(txt);
                        par.distortion_correction_distance = tmp;
                        txt = sprintf( '\nENTER OUTER OFFSET:');
                        tmp = input(txt);
                        par.distortion_correction_outer_offset = tmp;
                    case 'r'
                        %%
                        if exist('vol','var')
                            fprintf('SELECT ROI FOR METRIC CALCULATION:')
                            fprintf('\n 1. Set ROI using the mouse'  )
                            fprintf('\n 2. Double click on rectangle OR right click on rectangle and choose ''Crop Image'''  )
                            fprintf('\n 3. Close figure window'  )
                            pause(1)
                            figure('Name','Rotation axis offset: Set ROI for metric calculation')
                            z = round(size(vol,3) / 2);
                            im = vol(:,:,z);
                            im = FilterHisto(im,4);
                            im = normat(im);
                            [~,rect] = imcrop(im);
                            rect = round(rect);
                            fprintf('\n ROI selection: y0, x0, h, w = %u %u %u %u', rect)
                            x = rect(1) + rect(3)/2 - size(im,2)/2;
                            y = rect(2) + rect(4)/2 - size(im,1)/2;
                            d = sqrt(x^2 + y^2);
                            fprintf('\n ROI distance from center: %u pixel',round(d))
                            itomo.rot_axis_offset_metric_roi = rect;
                        else
                            fprintf('\nNo reco available. Run loop first to reconstruct a volume.')
                        end
                        
                        %%
                end % switch inp
                % Query parameters text
                txt = [...
                    '\n\nENTER RANGE OF ROTATION AXIS OFFSETS'...
                    '\n if empty: use default range, '...
                    '\n if scalar: use value & end interactive mode, '...
                    '\n if ''s'': change slice, '...
                    '\n if ''d'': enter debug mode, '...
                    '\n if ''c'': set distortion correction, '...
                    '\n if ''r'': set ROI for metric calculation, '...
                    '\n: '];
                inp = input( txt );
            end % while ischar( inp )
            if isempty( inp )
                offset = itomo.rot_axis_offset + interactive_mode.rot_axis_pos_default_search_range;
                fprintf( 'using default range\n' )
            else
                offset = inp;
            end  % isempty( inp )
            
            % Set offset or loop over offset range
            if isscalar( offset )
                fprintf( ' new rotation axis offset : %.2f (before: %.2f)', offset, itomo.rot_axis_offset)
                itomo.rot_axis_offset = offset;
                
                % Loop over tilt again?
                if ~isempty( tilt )
                    inp = '';
                    while isempty( inp )
                        inp = input( '\n\nENTER TILT LOOP AGAIN? (''y''/1,''n''/0) ');
                    end
                    switch lower( inp )
                        case {'y', 'yes', 1}
                            tilt = [];
                        case {'n', 'no', 0}
                    end
                end %  if ~isempty( tilt )
                
                % Loop over angles again?
                if ~isempty( angle_scaling )
                    inp = '';
                    while isempty( inp )
                        inp = input( '\n\nENTER ANGLE SCALING LOOP AGAIN? (''y''/1,''n''/0) ');
                    end
                    switch lower( inp )
                        case {'y', 'yes', 1}
                            angle_scaling = [];
                        case {'n', 'no', 0}
                    end
                end %  if ~isempty( tilt )
                
            else % when a range is given
                
                % Reco parameter
                [itomo.vol_shape, itomo.vol_size] = volshape_volsize( proj, itomo.vol_shape, itomo.vol_size, median(offset), 0);
                itomo.offset = offset;
                tomo.interactive_offset_range = offset;
                if ~isempty( tilt )
                    itomo.tilt = tilt;
                end
                if isscalar( angle_scaling )
                    itomo.angles = angle_scaling * angles + tomo.rot_angle_offset;
                else
                    itomo.angles = angles + tomo.rot_angle_offset;
                end
                
                % Reco
                [vol, metrics_offset] = find_rot_axis_offset(itomo,proj,par);
                
                % Metric minima
                [~, min_pos] = min(cell2mat({metrics_offset(:).val}));
                [~, max_pos] = max(cell2mat({metrics_offset(:).val}));
                
                % Print image number, rotation axis values, and different metrics
                fprintf( '\n no.' )
                fprintf( '%12s ', 'offset', metrics_offset.name)
                for nn = 1:numel(offset)
                    if offset(nn) == tomo.rot_axis_offset
                        fprintf('\n[\b%4u]\b%12.3f ', nn, offset(nn))
                    else
                        fprintf( '\n%4u%12.3f ', nn, offset(nn))
                    end
                    for mm = 1:numel(metrics_offset)
                        if min_pos(mm) == nn
                            fprintf( '[\b%12.2g]\b ', metrics_offset(mm).val(nn) )
                        elseif max_pos(mm) == nn
                            fprintf(2, '%12.2g ', metrics_offset(mm).val(nn) )
                        else
                            fprintf( '%12.2g ', metrics_offset(mm).val(nn) )
                        end
                    end
                end
                
                % Plot metrics
                h_rot_off = figure('Name', 'OFFSET: metrics', 'WindowState', window_state);
                ind = 2:7; %[1:4 6:7];
                Y = cell2mat({metrics_offset(ind).val});
                plot( offset, Y, '-+');
                ax1 = gca;
                axis tight
                xlabel( 'offset' )
                legend( metrics_offset(ind).name )
                ax2 = axes( 'Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
                x2 = 1:numel( offset );
                line( x2, 1, 'Parent', ax2 )
                ax2.XLim = x2([1 end]);
                %ax2.XTick = x2;
                ax2.XGrid = 'on';
                xlabel( 'index (image no.)' )
                set( ax1, 'YTick', [] ) % 'XTickMode', 'auto', 'XMinorTick', 'on')
                set( ax2, 'YTick', [] )
                title(sprintf('rotation axis: metrics VS offset'))
                drawnow
                saveas( h_rot_off, sprintf( '%s%s.png', write.fig_path, regexprep( h_rot_off.Name, '\ |:', '_') ) );
                
                % Display
                if show_stack_imagej
                    p = [write.interactive_path 'rot_axis_pos' filesep datestr(now, 'yyyymmddTHHMMSS') filesep];
                    mkdir(p)
                    if Bytes(vol(:,:,1)) < 2^32 - 1
                        fprintf('\nSaving tif image sequences:\n %s', p)
                        parfor nn = 1:size( vol, 3 )
                            filename = sprintf('%srot_axis_pos_%06u.tif', p, nn );
                            write32bitTIF(filename, vol(:,:,nn) );
                        end
                    else
                        fprintf('\nSaving h5 volume stack:\n %s', p)
                        filename = sprintf('%srot_axis_pos.h5', p );
                        if ~exist(filename,'file')
                            h5create(filename,'/volume',size(vol),'Datatype','single')
                        end
                        h5write(filename,'/volume',vol)
                        h5writeatt(filename,'/volume/','element_size_um', [1,1,1]);
                    end
                    p0 = pwd;
                    cd(p)
                    d = dir([p filesep '*.tif']);
                    while numel(d) ~= size(vol,3)
                        fprintf('\n WAITING FOR IMAGES TO BE APPEAR')
                        pause(1)
                        d = dir([p filesep '*.tif']);
                    end
                    if show_stack_imagej_use_virtual
                        if Bytes(vol(:,:,1)) < 2^32 - 1
                            fprintf('\nLoading tif sequence as virtual stack in imagej')
                            %unix('/asap3/petra3/gpfs/common/p05/jm/bin/imagej_opensequence &');
                            unix('imagej_opensequence &');
                        else
                            fprintf('\nLoading h5 volume as virtual stack in fiji')
                            %unix('/asap3/petra3/gpfs/common/p05/jm/bin/fiji_open_h5 &');
                            unix('fiji_open_h5 &');
                        end
                    else
                        if Bytes(vol(:,:,1)) < 2^32 - 1
                            fprintf('\nLoading all images in imagej (no virtual stack)')
                            %unix('/asap3/petra3/gpfs/common/p05/jm/bin/imagej_opensequence_novirt &');
                            unix('imagej_opensequence_novirt &');
                        else
                            fprintf('\nLoading h5 volume as virtual stack in fiji')
                            %unix('/asap3/petra3/gpfs/common/p05/jm/bin/fiji_open_h5_novirt &');
                            unix('fiji_open_h5_novirt &');
                        end
                    end
                    pause(2)
                    cd(p0)
                else
                    nimplay(vol, 1, [], 'OFFSET: sequence of reconstructed slices using different rotation axis offsets')
                end
            end % isscalar( offset )
            
            if isscalar( offset )
                
                %% TILT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if interactive_mode.rot_axis_tilt
                    
                    cprintf('RED', '\n\nEntering tilt loop:' )
                    cprintf('Magenta', '\n current rotation axis TILT : %g rad = %g deg', itomo.tilt, itomo.tilt * 180 / pi)
                    fprintf('\n calcul. rotation axis TILT : %g rad = %g deg', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
                    fprintf('\n image registration first and pi-shifted [index angle]: [%u %f] [%u %f] ',ind1,angles(ind1),ind2,angles(ind2))
                    fprintf('\n default tilt range : current TILT + [')
                    fprintf(' %.2g', interactive_mode.rot_axis_tilt_default_search_range )
                    fprintf(']' )
                    
                    % Loop over tilts
                    while ~isscalar( tilt )
                        
                        % Query parameters
                        inp = '';
                        while ischar( inp )
                            switch inp
                                case 's'
                                    slice = input( sprintf( '\nENTER ABSOLUTE [1,%u] OR RELATIVE [0,1] SLICE NUMBER : ', size( proj, 2 )) );
                                    if slice <= 1 && slice >= 0
                                        slice = round((size( proj, 2 ) - 1) * slice + 1 );
                                    end
                                    fprintf( ' new slice : %u (before: %u)', slice, itomo.slice );
                                    itomo.slice = slice;
                                case 'd'
                                    cprintf( 'RED', 'ENTERING DEBUG MODE. Continue with F5.' )
                                    keyboard
                            end % switch inp
                            
                            % Query parameters text
                            txt = [...
                                '\n\nENTER RANGE OF ROTATION AXIS TILTS'...
                                '\n if empty: use default range, '...
                                '\n if scalar: use value & end interactive mode, '...
                                '\n if ''s'': change slice, '...
                                '\n if ''d'': enter debug mode, '...
                                '\n: '];
                            inp = input( txt );
                        end % while ischar( inp )
                        if isempty( inp )
                            tilt = itomo.tilt + interactive_mode.rot_axis_tilt_default_search_range;
                            fprintf( 'using default range\n' )
                        else
                            tilt = inp;
                        end % isempty( inp )
                        
                        % Set tilt or loop over tilts
                        if isscalar( tilt )
                            fprintf( ' new rotation axis tilt : %f', tilt )
                            itomo.tilt = tilt;
                            
                        else
                            % Reco
                            itomo.tilt = tilt;
                            itomo.offset = offset;
                            if isscalar( angle_scaling )
                                itomo.angles = angle_scaling * angles + tomo.rot_angle_offset;
                            end
                            [vol, metrics_tilt] = find_rot_axis_tilt( itomo, proj);
                            
                            % Metric minima
                            [~, min_pos] = min( cell2mat( {metrics_tilt(:).val} ) );
                            [~, max_pos] = max( cell2mat( {metrics_tilt(:).val} ) );
                            
                            % Print image number and rotation axis tilt
                            fprintf( '\n no.' )
                            fprintf( '%12s ', 'tilt/rad', 'tilt/deg', metrics_tilt.name )
                            for nn = 1:numel(tilt)
                                fprintf( '\n%4u%12g %12g ', nn, tilt(nn), tilt(nn)/pi*180 )
                                for mm = 1:numel(metrics_tilt)
                                    if min_pos(mm) == nn
                                        fprintf( '[\b%12.3g]\b ', metrics_tilt(mm).val(nn) )
                                    elseif max_pos(mm) == nn
                                        fprintf( 2, '%12.3g ', metrics_tilt(mm).val(nn) )
                                    else
                                        fprintf( '%12.3g ', metrics_tilt(mm).val(nn) )
                                    end
                                end
                            end
                            
                            % Plot metrics
                            h_rot_tilt = figure('Name', 'TILT: metrics', 'WindowState', window_state);
                            ind = 6:7;
                            Y = cell2mat({metrics_tilt(ind).val});
                            plot( tilt, Y, '-+');
                            ax1 = gca;
                            axis tight
                            xlabel( 'tilt angle' )
                            legend( metrics_tilt(ind).name)
                            ax2 = axes( 'Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
                            x2 = 1:numel( tilt );
                            line( x2, 1, 'Parent', ax2 )
                            ax2.XLim = x2([1 end]);
                            %ax2.XTick = x2;
                            ax2.XGrid = 'on';
                            xlabel( 'index (image no.)' )
                            set( ax1, 'YTick', [] )
                            set( ax2, 'YTick', [] )
                            title(sprintf('rotation axis: metrics VS tilt'))
                            drawnow
                            saveas( h_rot_tilt, sprintf( '%s%s.png', write.fig_path, regexprep( h_rot_tilt.Name, '\ |:', '_') ) );
                            
                            % Display
                            if show_stack_imagej && Bytes(vol(:,:,1)) < 2^32 - 1
                                p = [write.interactive_path 'rot_axis_tilt' filesep datestr(now, 'yyyymmddTHHMMSS') filesep];
                                mkdir(p)
                                parfor nn = 1:size( vol, 3 )
                                    filename = sprintf('%srot_axis_pos_%06u.tif', p, nn );
                                    write32bitTIF(filename, vol(:,:,nn) );
                                end
                                p0 = pwd;
                                cd(p)
                                fprintf('\nLoading imagej')
                                if show_stack_imagej_use_virtual
                                    unix('/asap3/petra3/gpfs/common/p05/jm/bin/imagej_opensequence &');
                                    
                                else
                                    unix('/asap3/petra3/gpfs/common/p05/jm/bin/imagej_opensequence_novirt &');
                                end
                                pause(3)
                                cd(p0)
                            else
                                nimplay(vol, 1, [], 'TILT: sequence of reconstructed slices using different rotation axis tilts')
                            end
                        end
                        
                        if isscalar( tilt )
                            cprintf( 'RED', '\n\nCross check tilt by correlation of projections 0 and pi rad: ' )
                            fprintf( 'Registration of projection' )
                            tomo.rot_axis_position = im_shape_cropbin1 / 2 + tomo.rot_axis_offset;                            
                            % Compare projection at 0 pi and projection at 1 pi corrected for rotation axis tilt
                            y = round(size(proj,2)*[0.1 0.9]);
                            y = y(1):y(2);
                            s1 = size(proj,1);
                            if offset > 0
                                rpos = s1/2 + offset;
                                d = s1 - rpos;
                                pos = s1 - 2*d;
                                im1c = proj(pos:end,y,ind1);
                                im2c = flipud(proj(pos:end,y,ind2));
                            else
                                pos = round(2*(s1/2 + offset));
                                im1c = proj(1:pos,y,ind1);
                                im2c = flipud(proj(1:pos,y,ind2));
%                             im1c = RotAxisSymmetricCropping( proj(:,tomo.rot_axis_corr_area2,ind1),pos,1);
%                             im2c = flipud(RotAxisSymmetricCropping( proj(:,tomo.rot_axis_corr_area2,ind2),pos,1));
                            end
                            p = [write.im_path filesep 'correlation'];
                            CheckAndMakePath(p)                            
                            write32bitTIFfromSingle(sprintf('%s/proj%06u.tif',p,ind1),im1c)
                            write32bitTIFfromSingle(sprintf('%s/proj%06u.tif',p,ind2),im2c)
                            [optimizer, metric] = imregconfig('monomodal');
                            tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);
                            rot_axis_tilt_calc = asin( tform_calc.T(1,2) ) / 2;
                            %im2c_warped_calc =  imwarp(im2c, tform_calc, 'OutputView', imref2d(size(im1c)));
                            tform_int = tform_calc;
                            tform_int.T = [cos( 2 * itomo.tilt ) sin( 2 * itomo.tilt ) 0; ...
                                -sin( 2 * itomo.tilt ) cos( 2 * itomo.tilt ) 0 ; ...
                                0 0 1];
                                %tform_calc.T(3,1) tform_calc.T(3,2) 1];
                            % Remove translation if very large which is
                            % likely to be incorrect
%                             if abs(tform_int.T(3,1)) > 10
%                                 tform_int.T(3,1) = 0;
%                             end
%                             if abs(tform_int.T(3,2)) > 10
%                                 tform_int.T(3,2) = 0;
%                             end
                            fprintf( '\n registration matrix translation: x, y = %f %f', tform_int.T(3,1), tform_int.T(3,2) );
                            fprintf( '\n registration matrix scale: x, y = %f %f', tform_int.T(1,1), tform_int.T(2,2) );
                            fprintf( '\n registration matrix shear: x, y = %f %f', tform_int.T(1,2), tform_int.T(2,1) );
                            % Translation: [1 0 0, 0 1 0, t_x t_y 0]
                            % Scale: [s_x 0 0, 0 s_y 0, 0 0 1]
                            % Shear: [1 sh_y 0, sh_x 1 0, 0 0 1]
                            % Rotation: [cos(t) sin(t) 0, -sin(t) cos(t) 0, 0 0 1]
                            
                            Rin = imref2d(size(im1c));
                            Rin.XWorldLimits = Rin.XWorldLimits-1*mean(Rin.XWorldLimits);
                            Rin.YWorldLimits = Rin.YWorldLimits-1*mean(Rin.YWorldLimits);
                            %im2c_warped_int =  imwarp(im2c,tform_int,'OutputView',Rin);
                            %im2c_warped_int =  imwarp(im2c,Rin,tform_int);
                            fprintf( '\n Rotating 2nd image')
                            im2c_warped_int = imrotate(im2c,-2 * tilt * 180 / pi,'bilinear','crop');
                            %xt = 4 *ceil( 3 * abs( sin(2*itomo.tilt) ) * max( size(im1c)) ) + 10 + 2 *ceil(max( abs(tform_int.T(3,1:2))));
                            %if xt < size( im1c,1)  -10 && xt < size( im1c,2)  -10
                            fprintf( '\n current rotation axis tilt from interactive mode: %g rad (%g deg)', itomo.tilt, itomo.tilt * 180 / pi)
                            fprintf( '\n calcul. rotation axis tilt from registration    : %g rad (%g deg)', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)                            
                            x0 = max([1 + 2 * ceil(sin(2 * tilt) .* size(im1c,1)),1]);
                            y0 = max([1 + 2 * ceil(sin(2 * tilt) .* size(im1c,2)),1]);
                            x1 = size(im1c,1)-x0;
                            y1 = size(im1c,2)-y0;
                            x = x0:x1;
                            y = y0:y1;
                            im1 = im1c(x,y)';
                            im2org = im2c(x,y)';
                            im2 = im2c_warped_int(x,y)';                            
                            out = ImageCorrelation(im1,im2);
                            fprintf( '\n relative shift after rotation: %.1f, %0.1f',out.shift1,out.shift2)
                            fprintf( '\n Interpolating 2nd image on shifted grid')
                            [X,Y] = meshgrid(1:size(im1,2),1:size(im1,1));
                            im2i = interp2(X,Y,im2,X - out.shift2,Y - out.shift1,'linear',mean2(im2));                            
                            fprintf( '\n Cropping on common grid')
                            if out.shift2 > 0
                                y0 = 1 + ceil(abs(out.shift2));
                                y1 = size(im1,2);
                            else
                                y0 = 1;
                                y1 = size(im1,2) + floor(out.shift2);
                            end
                            if out.shift1 > 0
                                x0 = 1 + ceil(abs(out.shift1));
                                x1 = size(im1,1);
                            else
                                x0 = 1;
                                x1 = size(im1,1) + floor(out.shift1);
                            end
                            % tilt crop, worst case
                            d = ceil(max(size(im1)) * sin(abs(tilt))) + 1;
                            x = x0+d:x1-d;
                            y = y0+d:y1-d;
                            im1m = im1(x,y);
                            im2orgm = im2org(x,y);
                            im2m = im2(x,y);
                            im2im = im2i(x,y);
                            p2 = [write.im_path filesep 'correlation' filesep 'cropped'];
                            CheckAndMakePath(p2)   
                            write32bitTIFfromSingle(sprintf('%s/proj%06u.tif',p2,ind1),im1m)
                            write32bitTIFfromSingle(sprintf('%s/proj%06u.tif',p2,ind2),im2orgm)
                            write32bitTIFfromSingle(sprintf('%s/proj%06u_warped.tif',p2,ind2),im2m)
                            write32bitTIFfromSingle(sprintf('%s/proj%06u_warped_shifted.tif',p2,ind2),im2im)
                            name = sprintf( 'TILT: registered projections at %g and %g degree. rot axis tilt from INTERACTIVE mode: %g, rot axis offset: %g', angles(ind1)/pi*180, angles(ind2)/pi*180, tilt, offset);
                            nimplay(cat(3,im1m,im2im), 1, 0, name)
                            %nimplay( cat(3, im1c(xt:end-xt,xt:end-xt)', im2c_warped_int(xt:end-xt,xt:end-xt)'), 1, 0, name)
                            
%                             name = sprintf( 'TILT: registered projections at %g and %g degree. corrected. rot axis tilt from REGISTRATION: %g, rot axis offset: %g', angles(ind1)/pi*180, angles(ind2)/pi*180, rot_axis_tilt_calc, offset);
%                             %nimplay( cat(3, im1c(xt:end-xt,xt:end-xt)', im2c_warped_calc(xt:end-xt,xt:end-xt)'), 1, 0, name)
%                             nimplay( cat(3, im1c', im2c_warped_calc'), 1, 0, name)
%                             
                            inp = input( '\n\nENTER ROTATION AXIS TILT, if empty use current tilt: ');
                            if ~isempty( inp )
                                tilt = inp;
                            end
                            %end
                            
                            % Loop over offsets again?
                            inp = [];
                            while isempty( inp )
                                inp = input( '\nENTER OFFSET LOOP AGAIN? (''y''/1,''n''/0) ');
                            end
                            switch lower( inp )
                                case {'y', 'yes', 1}
                                    offset = [];
                                case {'n', 'no', 0}
                            end
                        end % if isscalar( tilt )
                    end % while ~isscalar( tilt )
                    itomo.tilt = tilt;
                end % if interactive_mode.rot_axis_tilt
                
                %% ANGLES
                if interactive_mode.angles
                    
                    cprintf( 'RED', '\n\nEntering angle loop:' )
                    ang_min = min( angles );
                    ang_max = max( angles );
                    ang_stride = mean( angles(2:end) - angles(1:end-1)  );
                    rel_ang_stride = ang_stride / ang_max;
                    if isempty( interactive_mode.angle_scaling_default_search_range )
                        r = 1 + rel_ang_stride*(-5:5);
                        interactive_mode.angle_scaling_default_search_range = r;
                    end
                    fprintf( '\n\n default angles : [min, max, average stride] = [%g %g %g]', ang_min, ang_max, ang_stride  )
                    fprintf( '\n default stride / max : %g', rel_ang_stride )
                    fprintf( '\n default scaling range : ' )
                    fprintf( ' %g', r )
                    
                    % Loop over angle_scalings
                    while ~isscalar( angle_scaling )
                        
                        % Query parameters
                        inp = '';
                        while ischar( inp )
                            switch inp
                                case 's'
                                    slice = input( sprintf( '\nENTER ABSOLUTE [1,%u] OR RELATIVE [0,1] SLICE NUMBER : ', size( proj, 2 )) );
                                    if slice <= 1 && slice >= 0
                                        slice = round((size( proj, 2 ) - 1) * slice + 1 );
                                    end
                                    fprintf( ' new slice : %u (before: %u)', slice, itomo.slice );
                                    itomo.slice = slice;
                                case 'd'
                                    cprintf( 'RED', 'ENTERING DEBUG MODE. Continue with F5.' )
                                    keyboard
                            end % switch inp
                            
                            % Query parameters text
                            txt = [...
                                '\n\nENTER (RANGE OF) ANGLE SCALING(S)'...
                                '\n if empty: use default range, '...
                                '\n if scalar: use value & end interactive mode, '...
                                '\n if ''s'': change slice, '...
                                '\n if ''d'': enter debug mode, '...
                                '\n: '];
                            inp = input( txt );
                        end % while ischar( inp )
                        if isempty( inp )
                            angle_scaling = interactive_mode.angle_scaling_default_search_range;
                            fprintf( 'using default range\n' )
                        else
                            angle_scaling = inp;
                        end % isempty( inp )
                        
                        % Set angles
                        if isscalar( angle_scaling )
                            fprintf( ' new angle scaling : %.2f (before: %.2f)', angle_scaling, tomo.angle_scaling )
                            itomo.angle_scaling = angle_scaling;
                            
                        else
                            % Reco
                            itomo.offset = offset;
                            itomo.tilt = tilt;
                            itomo.angle_scaling = angle_scaling;
                            [vol, metrics_angle_scaling] = find_angle_scaling( itomo, proj, angles );
                            
                            % Metric minima
                            [~, min_pos] = min( cell2mat( {metrics_angle_scaling(:).val} ) );
                            [~, max_pos] = max( cell2mat( {metrics_angle_scaling(:).val} ) );
                            
                            % Print image number and angle_scaling
                            fprintf( '\n no.' )
                            fprintf( '%12s ', 'scaling', metrics_angle_scaling.name )
                            for nn = 1:numel(angle_scaling)
                                fprintf( sprintf( '\n%4u%12f ', nn, angle_scaling(nn) ) )
                                for mm = 1:numel(metrics_angle_scaling)
                                    if min_pos(mm) == nn
                                        fprintf( '[\b%12.3g]\b ', metrics_angle_scaling(mm).val(nn) )
                                    elseif max_pos(mm) == nn
                                        fprintf( 2, '%12.3g ', metrics_angle_scaling(mm).val(nn) )
                                    else
                                        fprintf( '%12.3g ', metrics_angle_scaling(mm).val(nn) )
                                    end
                                end
                            end
                            
                            % Plot metrics
                            h_rot_angle_scaling = figure('Name', 'ANGLES: metrics', 'WindowState', window_state);
                            ind = 1:7;
                            Y = cell2mat({metrics_angle_scaling(ind).val});
                            plot( angle_scaling, Y, '-+');
                            ax1 = gca;
                            axis tight
                            xlabel( 'angle scaling' )
                            legend( metrics_angle_scaling(ind).name)
                            ax2 = axes( 'Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
                            x2 = 1:numel( angle_scaling );
                            line( x2, 1, 'Parent', ax2 )
                            ax2.XLim = x2([1 end]);
                            %ax2.XTick = x2;
                            ax2.XGrid = 'on';
                            xlabel( 'index (image no.)' )
                            set( ax1, 'YTick', [] )
                            set( ax2, 'YTick', [] )
                            title(sprintf('angles: metrics VS angle_scaling'))
                            drawnow
                            saveas( h_rot_angle_scaling, sprintf( '%s%s.png', write.fig_path, regexprep( h_rot_angle_scaling.Name, '\ |:', '_') ) );
                            
                            % Display
                            if show_stack_imagej && Bytes(vol(:,:,1)) < 2^32 - 1
                                p = [write.interactive_path 'rot_axis_angles' filesep datestr(now, 'yyyymmddTHHMMSS') filesep];
                                mkdir(p)
                                parfor nn = 1:size( vol, 3 )
                                    filename = sprintf('%srot_axis_pos_%06u.tif', p, nn );
                                    write32bitTIF(filename, vol(:,:,nn) );
                                end
                                p0 = pwd;
                                cd(p)
                                fprintf('\nLoading imagej')
                                if show_stack_imagej_use_virtual
                                    unix('/asap3/petra3/gpfs/common/p05/jm/bin/imagej_opensequence &');
                                    
                                else
                                    unix('/asap3/petra3/gpfs/common/p05/jm/bin/imagej_opensequence_novirt &');
                                end
                                pause(3)
                                cd(p0)
                            else
                                nimplay(vol, 1, [], 'ANGLES: sequence of reconstructed slices using different angle scalings')
                            end
                            
                        end
                        
                        if isscalar( angle_scaling )
                            % Loop over offsets again?
                            inp = [];
                            while isempty( inp )
                                inp = input( '\nENTER OFFSET LOOP AGAIN? (''y''/1,''n''/0) ');
                            end
                            switch lower( inp )
                                case {'y', 'yes', 1}
                                    offset = [];
                                case {'n', 'no', 0}
                            end
                        end % if isscalar( angle_scaling )
                        
                    end % while ~isscalar( angle_scaling )
                    tomo.angle_scaling = angle_scaling;
                    %tomo.angles = angle_scaling * angles + tomo.rot_angle_offset;
                    
                end % if interactive_mode.angles
                
            end % if isscalar( offset )
            
        end % while ~isscalar( offset )
        
        tomo.rot_axis_offset = offset;
        tomo.rot_axis_position = im_shape_cropbin1 / 2 + tomo.rot_axis_offset;
        
        % Save last sequence
        if exist( 'vol', 'var' ) && (Bytes(vol)/4 < 2^32 - 1)
            filename = sprintf( '%srot_axis_sequence.gif', write.fig_path );
            write_gif( vol, filename )
        end
        
        if ~lamino
            tomo.rot_axis_tilt_camera = itomo.tilt;
        else
            tomo.rot_axis_tilt_lamino = itomo.tilt;
        end
        
        tint = toc - tint;
        cprintf( 'RED', '\nEND OF INTERACTIVE MODE\n' )
    end % if interactive_mode.rot_axis_pos
    
    if isscalar( angle_scaling )
        angles = angle_scaling * angles;
    end
    
    ro = tomo.rot_axis_offset( ceil( numel( tomo.rot_axis_offset ) / 2 ) );
    rp = tomo.rot_axis_position( ceil( numel( tomo.rot_axis_position ) / 2 ) );
    fprintf( '\n rotation axis offset: %.2f', ro );
    fprintf( '\n rotation axis position: %.2f', rp );
    fprintf( '\n rotation axis tilt camera: %g rad (%g deg)', tomo.rot_axis_tilt_camera, tomo.rot_axis_tilt_camera * 180 / pi)
    fprintf( '\n rotation axis tilt lamino: %g rad (%g deg)', tomo.rot_axis_tilt_lamino, tomo.rot_axis_tilt_lamino * 180 / pi)
    [tomo.vol_shape, tomo.vol_size] = volshape_volsize( proj, tomo.vol_shape, tomo.vol_size, ro, 1);
    
    %% Display 0/pi projection: original and registered with tilt
    if par.visual_output && interactive_mode.rot_axis_tilt && lamino
        figure('Name','TILT: Projections at 0 and pi cropped symmetrically to rotation center', 'WindowState', window_state);
        n = 2;
        m = 2;
        
        subplot(m, n, 1)
        imsc1( im1c )
        axis equal tight
        title(sprintf('proj at 0'))
        colorbar
        
        subplot(m, n, 2)
        imsc1( im2c )
        axis equal tight
        title(sprintf('proj at pi'))
        colorbar
        
        xt = ceil( 3 * abs( sin(2*itomo.tilt) ) * max( size(im1c)) ) + 2;
        if xt > size( im1c,1)  -10 || xt > size( im1c,2)  -10
            xt = 1;
        end
        subplot(m, n, 3)
        imsc1( abs( im1c(xt:end-xt,xt:end-xt) - im2c(xt:end-xt,xt:end-xt) ) )
        axis equal tight
        title(sprintf('difference: original projections'))
        colorbar
        
        subplot(m, n, 4)
        [optimizer, metric] = imregconfig('monomodal');
        if ~exist( 'tform_int', 'var' )
            tform_int = imregtform(im2c, im1c, 'rigid', optimizer, metric);
        end
        im2c_warped_int =  imwarp(im2c, tform_int, 'OutputView',imref2d(size(im1c)));
        imsc1( abs( im1c(xt:end-xt,xt:end-xt) - im2c_warped_int(xt:end-xt,xt:end-xt) ) )
        axis equal tight
        title(sprintf('difference: registred projections'))
        colorbar
        
        drawnow
    end % if interactive_mode.rot_axis_tilt && par.visual_output
end % if tomo.run || tomo.run_interactive_mode
