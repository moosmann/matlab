function tomo = find_rot_axis_offset_auto(tomo, proj, par, write, interactive_mode)

offset = tomo.rot_axis_pos_search_range;
search_metric = tomo.rot_axis_pos_search_metric;

if ~isempty(offset) || ~isempty(search_metric)
        fprintf( '\nAutomatic rotation axis position determination.')
        fprintf( '\n metric : %s', search_metric )
        fprintf( '\n search range : ')
        fprintf( '%g ', offset )
        fprintf( '\n' )

    tomo_auto = tomo;
    tomo_auto.offset = offset;
    [~, metrics_offset] = find_rot_axis_offset( tomo_auto, proj);
    
    
    if tomo.rot_axis_pos_search_verbose
        
        % Metric minima
        [~, min_pos] = min(cell2mat({metrics_offset(:).val}));
        [~, max_pos] = max(cell2mat({metrics_offset(:).val}));
        
        % Print image number, rotation axis values, and different metrics
        fprintf( ' no.' )
        fprintf( '%11s', 'offset', metrics_offset.name)
        for nn = 1:numel(offset)
            
            cprintf( 'Black', '\n%4u%11.3f', nn, offset(nn))
            
            for mm = 1:numel(metrics_offset)
                if min_pos(mm) == nn
                    cprintf( 'Red', '%11.2g', metrics_offset(mm).val(nn) )
                elseif max_pos(mm) == nn
                    cprintf( 'Blue', '%11.2g', metrics_offset(mm).val(nn) )
                else
                    cprintf( 'Black', '%11.2g', metrics_offset(mm).val(nn) )
                end
            end
        end
    end
    
    % Plot metrics
    h_rot_off = figure('Name', 'OFFSET AUTO: metrics', 'WindowState', par.window_state);
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
    
    % pick metric
    m = strcmp( {metrics_offset.name}, search_metric );
    m = metrics_offset(m);
    m = m.val;
    % search for extrema and average
    switch tomo.rot_axis_pos_search_extrema
        case 'max' 
            v = max(m);
        case 'min'
            v = min(m);
    end    
    ind = (1:numel(m))';
    ind = ind( m == v );
    offset = mean( offset(ind) );
        
    fprintf( '\n offset found : %g', offset )
    tomo.rot_axis_offset = offset;
    
end

