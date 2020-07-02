% 30 oct 2014
ca
clear all

%% Figure parameter
% File and folder names
FigName = 'PhaseFilter_QP';
SavePath = sprintf('/mnt/tomoraid-LSDF/users/moosmann/phd/figures/%s/',FigName);
MakePath(SavePath)
% parameter
fontSize = 20;

%% Simulation parameter
edp = [EnergyConverter(1e-10); 0.1; 1e-6];
regpar = 1;
bf = 0.2;

%% Data
xdim = 4*2024;
pftie = PhaseFilter('tie',[xdim 1],edp,regpar,bf);
pfctf = PhaseFilter('ctf',[xdim 1],edp,regpar,bf);
pfqp = PhaseFilter('qp',[xdim 1],edp,regpar,bf);
pfqp2 = PhaseFilter('qp2',[xdim 1],edp,regpar,bf);
%pfqphs = PhaseFilter('qphalfsine',[xdim 1],edp,regpar,bf+3*a);

%% 2D Plots
% Create data

pfqp = sign(pfqp).*(abs(pfqp)-0.05);
xmax = floor(xdim/2);
xplot = 1:xmax;
X = xplot/xmax;
Y = [ pftie(xplot,1) pfctf(xplot,1)  pfqp(xplot,1) pfqp2(xplot,1) ];

% Create figure
figure1 = figure('Name','PhD thesis. Fourier space filte');

% Create axes
axes1 = axes('Parent',figure1,'FontSize',fontSize);
box(axes1,'on');
hold(axes1,'all');
set(axes1,'FontSize',fontSize)
% Create plot
plot1 = plot(X,Y,'Parent',axes1,'LineStyle','d','Marker','.');
% Plot line width
markSize = 8;
set(plot1(1),'Color','blue','MarkerSize',markSize);
set(plot1(2),'Color','red','MarkerSize',markSize);
set(plot1(3),'Color','green','MarkerSize',markSize+1);
set(plot1(4),'Color','black','MarkerSize',markSize);

% Figure properties
xlabel('xlabel','FontSize',fontSize)
ylabel('ylabel','FontSize',fontSize)

% Create legend
legend(axes1,'show');
legend('location','SouthWest');
legend('boxoff');

set(axes1, 'Box', 'off');
set(axes1,'FontSize',fontSize)

%% Save figure
saveas(figure1,sprintf('%s%s',SavePath,'Fig1.eps'),'epsc2')
saveas(gcf,sprintf('%s%s',SavePath,'Fig1_gcf.eps'),'epsc2')


