% 30 oct 2014
ca
clear all

%% Figure parameter
% File and folder names
FigName = 'TransferFunction';
SavePath = sprintf('/mnt/tomoraid-LSDF/users/moosmann/phd/figures/%s/',FigName);
MakePath(SavePath)
% parameter
fontSize = 20;

%% Simulation parameter
edp = [EnergyConverter(1e-10); 0.1; 1e-6];
regpar = 10;
bf = 0;

%% Data
xdim = 4*2024;
tftie = TransferFunction('tie',[xdim 1],edp,regpar,bf);
tfctf = TransferFunction('ctf',[xdim 1],edp,regpar,bf);


%% 2D Plots
% Create data
xmax = floor(xdim/2);
xplot = 1:xmax;
X = xplot/xmax;
Y = [tftie(xplot,1) tfctf(xplot,1)];

% Create figure
figure1 = figure('Name','PhD thesis. Fourier space filte');

% Create axes
axes1 = axes('Parent',figure1,'FontSize',fontSize);
box(axes1,'on');
hold(axes1,'all');
set(axes1,'FontSize',fontSize)
% Create plot
plot1 = plot(X,Y,'Parent',axes1);
% Plot line width
set(plot1,'LineWidth',2);
set(plot1(1),'Color','blue');
set(plot1(2),'Color','red');
% 
ylim(axes1,[-2 4]);

%ylim(plot1(1),[0 2]);
% Figure properties
xlabel('xlabel','FontSize',fontSize)
ylabel('ylabel','FontSize',fontSize)

%legend2 = legend(arrayfun(@(val) sprintf('z = %g m',val),z,'UniformOutput',0));
% Create legend
legend(axes1,'show');
legend('location','SouthWest');
legend('boxoff');

set(axes1, 'Box', 'off');
set(axes1,'FontSize',fontSize)

%% Save figure
saveas(figure1,sprintf('%s%s',SavePath,'Fig1.eps'),'epsc2')
saveas(gcf,sprintf('%s%s',SavePath,'Fig1_gcf.eps'),'epsc2')


