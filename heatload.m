% radiation induced temperature change measurements
%
% The tube was center with a metal screw as a position marker. Two
% thermistors can be inserted into the tube with motorized stage. The
% monochromatic x-ray beam of 1.5% bandwidth has average energy 31keV with
% 1mm Al filter in the upstream beamline. X-ray beam size was defined by
% the aperture by the end of the beam line to 3.2x3.2mm^2. The head of the
% thermistors are in the middle of the view. During rotation stage spinning
% the thermisters was moved up by 4mm so they were out of the beam during
% the spinning. The measurements were done as following.
%
% 1. measure temperature at the beam center level by inserting two
% thermistors in the tube. beam was off. this is the temperature before
% beam illumination; 100 temperature measurements were done in 20 seconds.
% 
% 2. move the thermistors out of the beam; turn on the beam; wait for 5
% seconds and then rotate the tube with stage at speed 18deg/s in range
% [-9, 369]. this takes about 22 seconds; take 100 temperature measurements
% on the same time in 20seconds at out-beam position.   
%
% 3. turn off the beam; move the thermistors back to the middle of the beam
% position; take 100 temperature measurements in 20 seconds. this is the
% temperature after x-ray exposing.  
%
% idle for 507 seconds and repeat the measurements again. Totally repeat 15
% times. 
%
% temperature reading from the amplifier refreshed every 0.1sec.
%
% See the screenshot of the the scanning picture for more detail time
% estimation. 

ttotal = 3*60+43+2/60; % in min
%% Read text file
InputPath = '/mnt/tomoraid-LSDF/users/moosmann/Nature/suppl_fig_3/';
SavePath = InputPath;
FileName = 'rad_induced_temp2';
% Open file.
fid = fopen([InputPath FileName],'r');
% Scan text in the file. The second next line uses a syntax which skips
% the rest of line '%*[^\n]'.
c = textscan(fid,'Data written 12 December 2012 %f:%f:%f%*[^\n]\n%*[^\n]\n%*s : %f%*[^\n]\n%*s : %f%*[^\n]\n%*[^\n]\n');
clear InputPath FileName fid d
%% 
nn = numel(c{1}(:));
trel = datenum([zeros(nn,1) zeros(nn,1) zeros(nn,1) c{1} c{2} c{3}]);
trel = trel - trel(1);
trel = single(trel/trel(6))/60;
% Temperature in Degree Celcius 
%Tborder = c{4} - 273.15;Tcenter = c{5} - 273.15;
Tborder = c{4};
Tcenter = c{5};
%% Average
pts = nn/5;
trel = mean(reshape(trel,[5 pts]));
Tborder = mean(reshape(Tborder,[5 pts]));
Tcenter = mean(reshape(Tcenter,[5 pts]));
%% Plot
MarkerSize = 5;
plotAllPoints = 0;
plotSeparately = 0;
%% All measurement points
if plotAllPoints(1)
    x = 1:1:pts;
    figure('Name','Central temperature VS time')
    plot(trel(x),Tcenter(x),'.','MarkerSize',MarkerSize)
    figure('Name','peripheral temperature VS time')
    plot(trel(x),Tborder(x),'.','MarkerSize',MarkerSize)
end
t = trel;
Tb = Tborder;
Tc = Tcenter;
% Substract mean
%Tb = Tb - mean(Tb(:));Tc = Tc - mean(Tc(:));
trel = reshape(trel,[20 pts/20]);
Tborder = reshape(Tborder,[20 pts/20]);
Tcenter = reshape(Tcenter,[20 pts/20]);
% Time
trel1 = trel(:,1:3:end);trel1 = mean(trel1);
trel2 = trel(:,2:3:end);trel2 = mean(trel2);
trel3 = trel(:,3:3:end);trel3 = mean(trel3);
% Central temperature
Tcenter1 = Tcenter(:,1:3:end);Tcenter1 = mean(Tcenter1);
Tcenter2 = Tcenter(:,2:3:end);Tcenter2 = mean(Tcenter2);
Tcenter3 = Tcenter(:,3:3:end);Tcenter3 = mean(Tcenter3);
% peripheral temperature
Tborder1 = Tborder(:,1:3:end);Tborder1 = mean(Tborder1);
Tborder2 = Tborder(:,2:3:end);Tborder2 = mean(Tborder2);
Tborder3 = Tborder(:,3:3:end);Tborder3 = mean(Tborder3);
%% Seperate measurement into befoer, during, and after exposure
if plotSeparately(1)
    x = 1:1:pts/3/20;
    figure('Name','Central temperature before VS time')
    plot(trel1(x),Tcenter1(x),'.','MarkerSize',MarkerSize)
    figure('Name','Central temperature during VS time')
    plot(trel2(x),Tcenter2(x),'.','MarkerSize',MarkerSize)
    figure('Name','Central temperature after VS time')
    plot(trel3(x),Tcenter3(x),'.','MarkerSize',MarkerSize)
    figure('Name','peripheral temperature before VS time')
    plot(trel1(x),Tborder1(x),'.','MarkerSize',MarkerSize)
    figure('Name','peripheral temperature during VS time')
    plot(trel2(x),Tborder2(x),'.','MarkerSize',MarkerSize)
    figure('Name','peripheral temperature after VS time')
    plot(trel3(x),Tborder3(x),'.','MarkerSize',MarkerSize)
end
ca

%% Absolute Temp VS time
x=1:1440;
X1 = t(x);
YMatrix1 = [Tc(x); Tb(x)];
% Create figure
figure1 = figure('Position',[1 1 1000 300],...
    'Name','Absolute Temperature VS time');
% Create axes
axes1 = axes('Parent',figure1);%,'YTickLabel',{'25.6','25.8','26.0','26.2','26.4','26.6'});
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 ceil(X1(end)/10)*10]);
box(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'Marker','.','LineStyle','none');
set(plot1(1),'Color',[1 0 0],'DisplayName','central');
set(plot1(2),'Color',[0 0 1],'DisplayName','peripheral');
% Create xlabel
xlabel('t [min]');
% Create ylabel
ylabel('T [K]','Interpreter','tex');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','SouthEast');
% Save figure
saveas(gcf,[SavePath 'suppl_fig_AbsoluteTemperature.eps'],'epsc2')
set(gcf,'PaperPositionMode','Auto')
saveas(gcf,[SavePath 'suppl_fig_3a.eps'],'epsc2')
saveas(gcf,[SavePath 'suppl_fig_3a.tif'],'tif')

%% Absolute Temp VS time, ROI
x=1:60;
X1 = t(x)*60;
YMatrix1 = [Tc(x); Tb(x)];
% Create figure
figure1 = figure('Position',[1 1 500 300],...
    'Name','Temperature VS time/s');
% Create axes
axes2 = axes('Parent',figure1);%,'YTickLabel',{'25.6','25.7','25.8','25.9','26.0','26.1','26.2','26.3','26.4','26.5','26.6'});
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes2,[0 75]);
box(axes2,'on');
hold(axes2,'all');
% Create multiple lines using matrix input to plot
plot2 = plot(X1,YMatrix1,'Parent',axes2,'Marker','.','LineStyle','none');
set(plot2(1),'Color',[1 0 0],'DisplayName','central');
set(plot2(2),'Color',[0 0 1],'DisplayName','peripheral');
% Create xlabel
xlabel('t [s]');
% Create ylabel
ylabel('T [K]','Interpreter','tex');
% Create legend
legend2 = legend(axes2,'show');
set(legend2,'Location','SouthEast');
% Save figure
saveas(gcf,[SavePath 'suppl_fig_AbsoluteTemperature_Section.eps'],'epsc2')
set(gcf,'PaperPositionMode','Auto')
saveas(gcf,[SavePath 'suppl_fig_3b.eps'],'epsc2')
saveas(gcf,[SavePath 'suppl_fig_3b.tif'],'tif')

% %% Absolute Temp VS measuring points
% x=1:1440;
% X1 = x;
% YMatrix1 = [Tc(x); Tb(x)];
% % Create figure
% figure2 = figure('Position',[1 1 1000 300],...
%     'Name','Temperature VS measurement points');
% % Create axes
% axes2 = axes('Parent',figure2);
% % Uncomment the following line to preserve the X-limits of the axes
% xlim(axes2,[0 1450]);
% box(axes2,'on');
% hold(axes2,'all');
% % Create multiple lines using matrix input to plot
% plot4 = plot(X1,YMatrix1,'Parent',axes2,'Marker','.','LineStyle','none');
% set(plot1(1),'Color',[1 0 0],'DisplayName','central');
% set(plot1(2),'Color',[0 0 1],'DisplayName','peripheral');
% % Create xlabel
% xlabel('measuring points');
% % Create ylabel
% ylabel('T [K]','Interpreter','tex');
% %ylabel('\Delta T  [^oC]');
% % Create legend
% legend1 = legend(axes2,'show');
% set(legend1,'Location','SouthEast');
% % Save figure
% saveas(gcf,[SavePath 'suppl_fig_AbsoluteTemperatureVSmeasuringPts.eps'],'epsc2')
% set(gcf,'PaperPositionMode','Auto')
% saveas(gcf,[SavePath 'suppl_fig_AbsoluteTemperatureVSmeasuringPts_resize.eps'],'epsc2')


%% Difference
x=59:60:1440;
X1 = t(x);
Y1 = Tc(x)-Tb(x);
% Create figure
figure3 = figure('Position',[1 1 500 300],...
    'Name','Temperature center - peripheral after exposure  VS measurement points');
% Create figure
axes3 = axes('Parent',figure3);%,'YTickLabel',{'0.15','0.16','0.17','0.18','0.19','0.20','0.21','0.22','0.23'});
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes3,[0 230]);
box(axes3,'on');
hold(axes3,'all');
% Create plot
plot3 = plot(X1,Y1,'Marker','.','LineStyle','none');
% Create xlabel
xlabel('t  [min]');
% Create ylabel
ylabel('\Delta T_{C-P}  [K]','Interpreter','tex');
% Save figure
saveas(gcf,[SavePath 'suppl_fig_TemperatureDifference.eps'],'epsc2')
set(gcf,'PaperPositionMode','Auto')
saveas(gcf,[SavePath 'suppl_fig_3c.eps'],'epsc2')
saveas(gcf,[SavePath 'suppl_fig_3c.tif'],'tif')

%% Other
% x=1:60;
% figure('Name','Temperature center - peripheral  VS time')
% plot(x,Tc(x)-Tb(x),'.')
% x=1:1440;
% figure('Name','Temperature center - peripheral  VS measurement points')
% plot(x,Tc(x)-Tb(x),'.')