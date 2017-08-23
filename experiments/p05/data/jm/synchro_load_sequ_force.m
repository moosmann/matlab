clear all;

input_path =  '/asap3/petra3/gpfs/p05/2017/data/11003950/raw';
font_size = 12;
outpath = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/spie';

% Read force values for radiography and tomography
for nn=19:-1:0
    filename = sprintf('%s/syn13_55L_Mg10Gd_12w_load_%02u/force.dat',input_path, nn);
    f{nn+1} = read_dat( filename );
end
Y = [f{:}];
ft = [f{1:2:end}];
fr = [f{2:2:end}];

% Create indices for radiography and tomography
nt = 1201;
nr = 301;
for nn = 10:-1:1
    xt{nn} = (1:1:nt) + (nn-1)*(nt+nr);
    xr{nn} = nt + (1:nr) + (nn-1)*(nt+nr);
end
xt = [xt{:}];
xr = [xr{:}];



% Create figure
fig = figure(1);

% Create axes
axes1 = axes('Parent',fig,'FontSize',font_size,'XMinorTick','off');
hold(axes1,'on');

% Create plot
%plot( Y, '-', 'LineWidth',3);
plot( xt, ft, 'b.', xr, fr, 'r.', 'LineWidth',3);

% Create xlabel
xlabel('number of acquired images', 'FontSize',font_size);

% Create ylabel
ylabel('Force / N', 'FontSize',font_size);

% Set the remaining axes properties
box(axes1,'on');
axis(axes1,'tight');
set(axes1,'FontSize',font_size);

name = 'plot_force';

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig, sprintf( '%s/%s.png', outpath, name), 'png')
print(fig, sprintf( '%s/%s.pdf', outpath, name), '-dpdf')

fig.PaperUnits = 'centimeters';
fig.PaperSize = [12 6];
saveas(fig, sprintf( '%s/%s_2.png', outpath, name), 'png')
print(fig, sprintf( '%s/%s_2.pdf', outpath, name), '-dpdf')

if 0
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition',[0 0 1 1] );
set(gcf, 'PaperPositionMode', 'auto');


% Save figure
saveas(gcf, sprintf( '%s/%s.png', outpath, name), 'png')

saveas(gcf, sprintf( '%s/%s.jpg', outpath, name), 'jpeg')
saveas(gcf, sprintf( '%s/%s.eps', outpath, name), 'epsc')
saveas(gcf, sprintf( '%s/%s.svg', outpath, name), 'svg')
saveas(gcf, sprintf( '%s/%s.pdf', outpath, name), 'pdf')
end