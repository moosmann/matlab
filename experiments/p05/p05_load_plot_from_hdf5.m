function p05_load_plot_from_hdf5( filename, out_path, adc2force)

% Read out wait_force values of 4D load tomography data.
%
% filename : string, full filename
% out_path : string, default:''. directory where figures will be saved. if
%   empty uses (processed folder)/scan_name (without trailing indices)
% adc2force : scalar, default 4.8. Conversion from voltage read out of the
%   load cell to Newton

%% DEFAULT ARGUMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    filename = '/asap3/petra3/gpfs/p05/2018/data/11005553/raw/test/easyform_15x2N_relax60min_nexus.h5';
end
if nargin < 2
    out_path = '';
end
if nargin < 3
    adc2force = 10.4367;
    fprintf( '\nLoad cell conversion from Volt to Newtone set per default to %.f N. Check value!', adc2force)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf( '\n Plot load cell values')

[folder, scan_name] = fileparts( filename );

if isempty( out_path )
    out_path = regexprep( folder,'raw','processed');
end
CheckAndMakePath( out_path );

font_size = 12;
fprintf( '\n filename : %s', filename)
fprintf( '\n output path : %s', out_path)

%% Read load sequence data
time = double( h5read( filename, '/entry/hardware/adc1/Value/time') )';
value = double( h5read( filename, '/entry/hardware/adc1/Value/value') )';
time = time(2:end-1);
value = adc2force * value(2:end-1);

% Write CSV file
csvfilename = regexprep( filename , {'.h5', 'raw'}, {'.csv', 'processed'} );
fprintf( '\n create csv file : %s', csvfilename )
cHeader = {'time' 'force'}; %dummy header
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

%write header to file
fid = fopen( csvfilename,'w');
fprintf( fid,'%s\n',textHeader);
fclose( fid );

time = time / 3600 / 1e3;

datvalues(:,1) = time - time(1);
datvalues(:,2) = value;

%write data to end of file
dlmwrite(csvfilename,datvalues,'-append');

clear datvalues ;

%% Show and Save figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create figure
name = regexprep( scan_name,'\_',' ');
fig_name = sprintf( 'Force value measurements: %s', name );
fig = figure( 'Name', fig_name );

% Create axes
axes1 = axes('Parent',fig,'FontSize',font_size,'XMinorTick','off');
%hold(axes1,'on');

%plot( wait_force, 'b.', 'LineWidth',12);
plot( time - time(1), value, 'r.', 'LineWidth',12);
%legend( {'Waiting'}, 'FontSize', font_size,'Location','northwest')

% Create ylabel
xlabel('time / h', 'FontSize',font_size);
ylabel('Force / N', 'FontSize',font_size);

title( sprintf( '%s', name ) )

axis tight

% Set the remaining axes properties
box(axes1,'on');
axis(axes1,'tight');
set(axes1,'FontSize',font_size);

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

filename_fig = sprintf( '%s/plot_force_%s.png', out_path, scan_name);
fprintf( '\n figure filename : %s', filename )
saveas(fig, filename_fig, 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nFINISHED\n')
