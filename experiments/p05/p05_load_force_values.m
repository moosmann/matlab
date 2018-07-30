function p05_load_force_values( raw_path, scan_name, steps, out_path )

% Read out force values of 4D load tomography data.
%
% raw_path : path to raw folder
% scan_name : string, partial name of the sequence that all steps have in
%   common i.e. without trailing indices, underscores, etc
% steps : scalar or 1D vector, default: [] uses all. steps which will be
%   processed. if scalar, stops at number 'steps'
% out_path : string, default:''. directory where figures will be saved. if
%   empty uses: (processed folder)/scan_name

%% DEFAULT ARGUMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    raw_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/raw';
end
if nargin < 2
    %scan_name = 'syn002_6L_PEEK_4w';
    %scan_name = 'syn003_92L_Mg10Gd_4w';
    %scan_name = 'syn004_84L_Mg10Gd_4w';
    %scan_name = 'syn005_81L_Mg5Gd_8w';
    %scan_name = 'syn006_75R_Mg10Gd_8w';
    %scan_name = 'syn007_94L_Mg10Gd_8w';
    %scan_name = 'syn008_76R_Mg10Gd_8w';
    %scan_name = 'syn009_32R_PEEK_8w';
    %scan_name = 'syn010_19R_PEEK_4w';
    scan_name = 'syn011_14R_PEEK_4w';
    %scan_name = 'syn012_79L_Mg5Gd_8w';
end
if nargin < 3
    steps = []; % # of tomos to process, use low number for testing
end
if nargin < 4
    out_path = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf( '\n LOAD SEQUENCE PROCESSING')

scan_path = [raw_path filesep scan_name];

if isempty( out_path )
    out_path = regexprep( scan_path,'\raw','processed');
end
CheckAndMakePath( out_path );

%% Check value
adc2force = 4.8;
font_size = 12;
fprintf( '\n scan path : %s', scan_path)
fprintf( '\n scan name : %s', scan_name)

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

p = [raw_path filesep struct_scans(3).name];

%% Read in load sequence
time = [];
t1 = [];
timet = [];
force = [];
forcet = [];
pusher = [];
for nn = 1:num_steps
    ss = steps(nn);
    fprintf( '\n %2u : %s.', nn, struct_scans(ss).name)
    p = [raw_path filesep struct_scans(ss).name];
    filename = sprintf( '%s/%stime_adc1_p05pusher_wait.dat',p, scan_name);
    %fprintf( '\n%s %u', filename, exist(filename, 'file'))
    
    % Waiting
    if exist( filename, 'file')
        fw = read_dat( filename );
        t1 = [t1 fw(1,end-1)];
        time = [time fw(1,2:end-1)];
        force = [force adc2force * fw(2,2:end-1)];
        pusher = [pusher fw(3,2:end-1)];
    end
    
    % tomogram
    filename = sprintf( '%s/%s_nexus.h5', p, struct_scans(ss).name );
    if exist( filename, 'file')
        %h5i = h5info( filename );
        t = double( h5read( filename, '/entry/hardware/adc1/Value/time') );
        v = double( h5read( filename, '/entry/hardware/adc1/Value/value') );
        timet = [timet t(2:end-1)'];
        forcet = [forcet adc2force * v(2:end-1)'];
    end  
end

time = time / 60 / 60;
timet = timet / 3600 / 1e3;
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

if ~isempty( time ) && ~isempty( timet )
    plot( time - time(1), force, 'b.', timet - time(1), forcet, 'r.', 'LineWidth',12);
    legend( {'Waiting', 'Tomogram'}, 'FontSize', font_size,'Location','northwest')
else
    if ~isempty( time )
        plot( time - time(1), force, 'b.', 'LineWidth',12);
        legend( {'Waiting'}, 'FontSize', font_size,'Location','northwest')
    else
        plot( timet - timet(1), forcet, 'r.', 'LineWidth',12);
        legend( {'Tomogram'}, 'FontSize', font_size,'Location','northwest')
        
    end
end

% Create ylabel
ylabel('Force / N', 'FontSize',font_size);

% Create xlabel
xlabel('time / h', 'FontSize',font_size);

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
fprintf('\n')
