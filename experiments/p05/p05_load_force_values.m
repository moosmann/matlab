function p05_load_force_values( raw_path, scan_name, steps, out_path, adc2force, readhdf5 )

% Read out wait_force values of 4D load tomography data.
%
% raw_path : path to raw folder
% scan_name : string, partial name of the sequence that all steps have in
%   common i.e. without trailing indices, underscores, etc
% steps : scalar or 1D vector, default: [] uses all. steps which will be
%   processed. if scalar, stops at number 'steps'
% out_path : string, default:''. directory where figures will be saved. if
%   empty uses (processed folder)/scan_name (without trailing indices)
% adc2force : scalar, default 4.8. Conversion from voltage read out of the
%   load cell to Newton

%% DEFAULT ARGUMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    % e.g.
    raw_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/raw';
end
if nargin < 2
    % e.g.
    scan_name = 'syn011_14R_PEEK_4w';
end
if nargin < 3
    steps = []; % # of tomos to process, use low number for testing
end
if nargin < 4
    out_path = '';
end
if nargin < 5
    adc2force = 4.8;
    fprintf( '\nLoad cell conversion from Volt to Newtone set per default to %.f N. Check value!', adc2force)
end
if nargin < 6
    readhdf5 = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf( '\n LOAD SEQUENCE PROCESSING')

scan_path = [raw_path filesep scan_name];

if isempty( out_path )
    out_path = regexprep( scan_path,'raw','processed');
end
CheckAndMakePath( out_path );

font_size = 12;
fprintf( '\n scan path : %s', scan_path)
fprintf( '\n scan name : %s', scan_name)
fprintf( '\n output path : %s', out_path)

% scans to process
struct_scans = dir([scan_path '*']);
if numel( struct_scans ) == 0
    error( 'No scans found matching folder pattern: \n%s\n', scan_path )
end

if isempty( steps )
    steps = 1:numel( struct_scans);
end
if isscalar( steps )
    steps = 1:steps;
end
num_steps = numel( steps );

%p = [raw_path filesep struct_scans(3).name];

%% Read in load sequence
wait_time = [];
wait_time1 = [];
tomo_time = [];
tomo_time1 = [];
wait_force = [];
wait_adc = []; % legacy support for first scans
tomo_force = [];
wait_pusher = [];
tomo_pusher = [];
setforce_time = [];
setforce_force = [];

for nn = 1:num_steps
    ss = steps(nn);
    scan_name_nn = struct_scans(ss).name;
    fprintf( '\n %2u : %s.', nn, scan_name_nn)
    p = [raw_path filesep scan_name_nn];
    
    filename = sprintf( '%s/%stime_adc1_p05pusher_wait.dat',p, scan_name);
    % Waiting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist( filename, 'file')
        fprintf( 'time_adc1_p05pusher_wait.dat' )
        fw = read_dat( filename );
        wait_time1 = [wait_time1 fw(1,end-1)];
        wait_time = [wait_time fw(1,2:end-1)];
        wait_force = [wait_force adc2force * fw(2,2:end-1)];
        wait_pusher = [wait_pusher fw(3,2:end-1)];
        
        % Write CVS file
        csvfilename=[out_path filesep struct_scans(nn).name '_waiting.csv'];
        cHeader = {'wait_time' 'wait_force' 'wait_pusher'}; %dummy header
        commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
        commaHeader = commaHeader(:)';
        textHeader = cell2mat(commaHeader); %cHeader in text with commas
        
        %write header to file
        fid = fopen(csvfilename,'w');
        fprintf(fid,'%s\n',textHeader);
        fclose(fid);
        
        datvalues(:,1)=(wait_time-wait_time(1)) / 60 / 60;
        datvalues(:,2)=wait_force;
        datvalues(:,3)=wait_pusher;
        
        %write data to end of file
        dlmwrite(csvfilename,datvalues,'-append');
        
        clear datvalues ;
    end
    
    filename = sprintf( '%s/force.dat',p);
    if exist( filename, 'file' )        
        
        fw = read_dat( filename );
        wait_force = [wait_force fw(2:end-1)];
        
        filename = sprintf( '%s/posadc1.dat',p);
        fw = read_dat( filename );
        wait_adc = [wait_adc fw(2:end-1)];
        
        filename = sprintf( '%s/pospusher.dat',p);
        fw = read_dat( filename );
        wait_pusher = [wait_pusher fw(2:end-1)];
        
        fprintf( ' force adc1 pusher dat' )
    end
    
    % tomogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename = sprintf( '%s/%stime_adc1_p05pusher_tomo.dat',p, scan_name_nn);
    if exist( filename, 'file')
        fprintf( 'time_adc1_p05pusher_tomo.dat' )
        fw = read_dat( filename );
        tomo_time1 = [tomo_time1 fw(1,end-1)];
        tomo_time = [tomo_time fw(1,2:end-1)];
        tomo_force = [tomo_force adc2force * fw(2,2:end-1)];
        tomo_pusher = [tomo_pusher fw(3,2:end-1)];
        
        % Write CVS file
        csvfilename=[out_path filesep struct_scans(nn).name '_tomo.csv'];
        cHeader = {'tomo_time' 'tomo_force' 'tomo_pusher'}; %dummy header
        commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
        commaHeader = commaHeader(:)';
        textHeader = cell2mat(commaHeader); %cHeader in text with commas
        
        %write header to file
        fid = fopen(csvfilename,'w');
        fprintf(fid,'%s\n',textHeader);
        fclose(fid);
        
        datvalues(:,1) = (tomo_time - tomo_time(1) ) / 60 / 60;
        datvalues(:,2) = tomo_force;
        datvalues(:,3) = tomo_pusher;
        
        %write data to end of file
        dlmwrite(csvfilename,datvalues,'-append');
        
        clear datvalues ;
    end
    
    filename = sprintf( '%s/%s_nexus.h5', p, struct_scans(ss).name );
    if readhdf5 && exist( filename, 'file')
        fprintf( ' nexus.h5' )
        %h5i = h5info( filename );
        t = double( h5read( filename, '/entry/hardware/adc1/Value/time') );
        v = double( h5read( filename, '/entry/hardware/adc1/Value/value') );
        tomo_time = [tomo_time t(2:end-1)'];
        tomo_force = [tomo_force adc2force * v(2:end-1)'];
        
        % Write CSV file
        csvfilename=[out_path filesep struct_scans(nn).name '_tomo.csv'];
        cHeader = {'wait_time' 'wait_force'}; %dummy header
        commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
        commaHeader = commaHeader(:)';
        textHeader = cell2mat(commaHeader); %cHeader in text with commas
        
        %write header to file
        fid = fopen(csvfilename,'w');
        fprintf(fid,'%s\n',textHeader);
        fclose(fid);
        
        datvalues(:,1)=(tomo_time-tomo_time(1)) / 3600 / 1e3;
        datvalues(:,2)=tomo_force;
        
        %write data to end of file
        dlmwrite(csvfilename,datvalues,'-append');
        
        clear datvalues ;
    end
    
    %% Set force and relaxation time
    s2 = dir( [p '_*'] );
    if numel(s2)
        p2 = [s2.folder filesep s2.name];
        filename = sprintf( '%s/%s_nexus.h5', p2, s2.name );
        
        if readhdf5 && exist( filename, 'file')
            fprintf( ' nexus.h5' )
            %h5i = h5info( filename );
            t = double( h5read( filename, '/entry/hardware/adc1/Value/time') );
            v = double( h5read( filename, '/entry/hardware/adc1/Value/value') );
            setforce_time = [setforce_time t(2:end-1)'];
            setforce_force = [setforce_force adc2force * v(2:end-1)'];
            
            % Write CSV file
            csvfilename=[out_path filesep struct_scans(nn).name '_tomo.csv'];
            cHeader = {'wait_time' 'wait_force'}; %dummy header
            commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
            commaHeader = commaHeader(:)';
            textHeader = cell2mat(commaHeader); %cHeader in text with commas
            
            %write header to file
            fid = fopen(csvfilename,'w');
            fprintf(fid,'%s\n',textHeader);
            fclose(fid);
            
            datvalues(:,1)=(tomo_time-tomo_time(1)) / 3600 / 1e3;
            datvalues(:,2)=tomo_force;
            
            %write data to end of file
            dlmwrite(csvfilename,datvalues,'-append');
            
            clear datvalues ;
        end
    end
end

wait_time = wait_time / 60 / 60;
setforce_time = setforce_time / 3600 / 1e3;
tomo_time = tomo_time / 3600 / 1e3;

%% Save
%filename = sprintf( '%s/%s_loadSequ_x%04u.gif', scan_path, scan_name, xx );
%fprintf( '\n output file: %s', filename)

%% Show and Save figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create figure
name = sprintf( 'Force value measurements: %s', scan_name );
fig = figure( 'Name', name );

% Create axes
axes1 = axes('Parent',fig,'FontSize',font_size,'XMinorTick','off');
%hold(axes1,'on');

if ~isempty( setforce_force )
    
    plot( tomo_time - tomo_time(1), tomo_force, 'r.', setforce_time - tomo_time(1), setforce_force, 'b.', 'LineWidth',12);
    legend( {'Tomogram', 'Set force / Waiting'}, 'FontSize', font_size,'Location','northwest')
else
    if isempty( wait_time ) && isempty( tomo_time )
        
        plot( wait_force, 'b.', 'LineWidth',12);
        legend( {'Waiting'}, 'FontSize', font_size,'Location','northwest')
        
        xlabel('projections', 'FontSize',font_size);
        
    else
        if isempty( wait_time ) || isempty( tomo_time )
            if ~isempty( wait_time )
                plot( wait_time - wait_time(1), wait_force, 'b.', 'LineWidth',12);
                legend( {'Waiting'}, 'FontSize', font_size,'Location','northwest')
            else
                plot( tomo_time - tomo_time(1), tomo_force, 'r.', 'LineWidth',12);
                legend( {'Tomogram'}, 'FontSize', font_size,'Location','northwest')
            end
        else
            plot( wait_time - tomo_time(1), wait_force, 'b.', tomo_time - tomo_time(1), tomo_force, 'r.', 'LineWidth',12);
            legend( {'Waiting', 'Tomogram'}, 'FontSize', font_size,'Location','northwest')
        end
        
        % Create xlabel
        xlabel('waiting / h', 'FontSize',font_size);
    end
end
% Create ylabel
ylabel('Force / N', 'FontSize',font_size);


plot_title = regexprep(scan_name,'\_',' ');
title( sprintf( '%s', plot_title ) )

axis tight

% Set the remaining axes properties
box(axes1,'on');
axis(axes1,'tight');
set(axes1,'FontSize',font_size);

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

filename = sprintf( '%s/plot_force_%s.png', out_path, scan_name);
saveas(fig, filename, 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nFINISHED\n')
