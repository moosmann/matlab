function tomo = find_rot_axis_offset_auto(tomo, proj, par, write, interactive_mode)
% Routine to automatically determine the rotation axis position by finding
% the extrema of a metric.

%% move variable check to main script
if tomo.rot_axis_search_auto == 0
    return
end

if isempty(tomo.rot_axis_search_range)
   tomo.rot_axis_search_range = interactive_mode.rot_axis_search_range;
end
if tomo.rot_axis_search_range_from_interactive
    tomo.rot_axis_search_range = tomo.interactive_offset_range;
end

offset = tomo.rot_axis_search_range;
search_metric = tomo.rot_axis_search_metric;
extrema = tomo.rot_axis_search_extrema;
while isempty(offset)
    cprintf( 'red', 'Rotation axis search range is undefined. Please enter search range:' );
    tmp = input( ' ' );
    if isvector( tmp ) && ~isscalar(tmp)
        offset = tmp;
        tomo.rot_axis_search_range = tmp;
    end
end
while isempty(search_metric)
    cprintf( 'red', 'Rotation axis search metric is undefinded. Enter metric (''neg'',''entropy'',''iso-grad'',''laplacian'',''entropy-ML'',''abs''):' );
    tmp = input( ' ' );
    if ischar( tmp )
        if sum(strcmp( tmp, {'neg','entropy','iso-grad','laplacian','entropy-ML','abs'} )) == 1
            search_metric = tmp;
            tomo.rot_axis_search_metric = tmp;
        end
    end
end
while isempty(extrema)
    cprintf( 'red', 'Rotation axis search extrema is undefinded. Enter extrema (''max'' or ''min''):' );
    tmp = input( ' ' );
    if ischar( tmp )
        if sum(strcmp( tmp, {'max','min'} )) == 1
            extrema = tmp;
           tomo.rot_axis_search_extrema = tmp;
        end
    end
end
if isempty(tomo.rot_axis_search_slice)
    tomo.rot_axis_search_slice = tomo.slice;
end
tomo.slice = tomo.rot_axis_search_slice;

if ~isempty(offset) || ~isempty(search_metric)
        fprintf( '\nAutomatic rotation axis position determination.')
        fprintf( '\n metric : %s', search_metric )
        fprintf( '\n search range : ')
        fprintf( '%g ', offset )
        fprintf( '\n' )

    tomo_auto = tomo;
    tomo_auto.offset = offset;
    [~, metrics_offset] = find_rot_axis_offset( tomo_auto, proj, par);
    
    if tomo.rot_axis_search_verbose
        % Metric minima
        [~, min_pos] = min(cell2mat({metrics_offset(:).val}));
        [~, max_pos] = max(cell2mat({metrics_offset(:).val}));
        % Print image number, rotation axis values, and different metrics
        fprintf( ' no.' )
        fprintf( '%11s', 'offset', metrics_offset.name)
        cp = 0;
        for nn = 1:numel(offset)
            if cp
                cprintf( 'Black', '\n%4u%11.3f', nn, offset(nn))
                cp = 0; pause(0.01)
            else
                fprintf( '\n%4u%11.3f', nn, offset(nn))
            end            
            for mm = 1:numel(metrics_offset)
                if min_pos(mm) == nn
                    cprintf( 'Red', '%11.2g', metrics_offset(mm).val(nn) )
                    cp = 1;
                elseif max_pos(mm) == nn
                    cprintf( 'Blue', '%11.2g', metrics_offset(mm).val(nn) )
                    cp = 1; pause(0.01)
                else
                    if cp
                        cprintf( 'Black', '%11.2g', metrics_offset(mm).val(nn) )
                        cp = 0; pause(0.01)
                    else
                        fprintf( '%11.2g', metrics_offset(mm).val(nn) )
                    end
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
    % search for extrema and average if there are multiple
    switch extrema
        case 'max' 
            v = max(m);
        case 'min'
            v = min(m);
    end
    ind = (1:numel(m))';
    ind = ind( m == v );
    metric_mean = mean(m(ind));
    offset_mean = mean( offset(ind) );

    
    % Fit metric
    if tomo.rot_axis_search_fit
        fprintf( '\n Fit metric %s with polynomial\n', search_metric )
        [f, gof] = fit(offset(:),m(:),'poly4');
        disp(gof)
        % find fit min
        foffset = min(offset):(max(offset)-min(offset))/numel(offset)/100:max(offset);
        switch extrema
            case 'max'
                [fit_val,fit_pos] = max( f(foffset) );
            case 'min'
                [fit_val,fit_pos] = min( f(foffset) );
        end
        fit_pos = foffset( fit_pos );
        h_fit = figure('Name', 'OFFSET AUTO: fit metric', 'WindowState', par.window_state);
        plot(f,offset,m);
        hold on
        plot( offset_mean, metric_mean, 'o', 'MarkerSize', 30)
        hold on
        plot( fit_pos, fit_val, 'x', 'MarkerSize', 30)
        set( gca, 'YTick', [] )
        axis tight
        xlabel( 'offset' )
        ylabel( 'metric' )
        legend( {search_metric, 'fit', sprintf('metric %s', extrema), sprintf('fit %s', extrema)} )
        title(sprintf('rotation axis: metric %s VS offset', search_metric))
        drawnow
        saveas( h_fit, sprintf( '%s%s.png', write.fig_path, regexprep( h_fit.Name, '\ |:', '_') ) );
        
        offset = fit_pos;
    else
        offset = offset_mean;
    end
            
    fprintf( '\n offset found : %g', offset )
    tomo.rot_axis_offset = offset;
    tomo.rot_axis_position = size(proj,1) / 2 + offset;
    
end

