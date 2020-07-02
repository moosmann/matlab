%Prior to evaluation the script 'CTFanlayisisXeno4cell.m' needs to be run or the
%data.
clear X1 X2 YMatrix1 plot axes1 axes2 figure1 figure2 plot1 plot2
X1=x;
SineArgPreFac = pi*EnergyConverter(energy)*distance/(paddim*pixelsize)^2;
figure1 = figure('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98404194812 29.67743169791],...
    'Name','Angular integrated line plot of FT of intensity contrast VS frequency');
% Create axes
X2 = SineArgPreFac*X1.^2;
%  axes1 = axes('Parent',figure1, ...
%       'YTickLabel',[0:400:1200],'FontSize',18);
%'XTick',[25 50 75 100 125 150 175], ...
axes1 = axes('Parent',figure1,'YTickLabel',{'','200','400','600','800','1000','1200','1400','1600','1800'},'FontSize',18);
% xlim(axes1,[25 180]);
% ylim(axes1,[0 1200]);
box(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
YMatrix1(4,:) = cf{2}(X1);
YMatrix1(3,:) = cf{1}(X1);
YMatrix1(2,:) = y(X1,2);
YMatrix1(1,:) = y(X1,1);
% Some options
set(gcf,'defaulttextinterpreter','none')
set(gcf,'position',[0 0, 1600 900]) 
% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',2);
set(plot1(1),'MarkerSize',15,'Marker','.','LineStyle','none','Color',[0 0 0],'LineWidth',0.5);
set(plot1(2),'MarkerSize',8,'Marker','x','LineStyle','none','Color',[0 0 1],'LineWidth',0.5);
set(plot1(3),'LineStyle','-','Color',[0 0 0]);
set(plot1(4),'LineStyle','--','Color',[0 0 1]);
% Create xlabel
xlabel('x','FontSize',22);
ylabel('y','FontSize',22);
% Save figure
saveas(gcf,'/home/moosmann/oe2011HMB/Fig-test_quadr.eps','epsc2')
set(gcf,'PaperPositionMode','Auto') 
saveas(gcf,'/home/moosmann/oe2011HMB/Fig-test.eps','epsc2')
saveas(gcf,'/home/moosmann/oe2011HMB/IntegratedLinePlotOfAbsIntensityOfXeno4cell.eps','epsc2')
% figure2 = figure('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%     'PaperSize',[20.98404194812 29.67743169791],...
%     'Name','Ratio of angular integrated line plot of FT of intensity contrast VS frequency');
% set(gcf,'position',[0 0, 1600 900]) 
% % Create multiple lines using matrix input to plot
% axes2 = axes('Parent',figure2,'FontSize',18);
% box(axes2,'on');
% hold(axes2,'all');
% plot2 = plot(X1,YMatrix1(2,:)./YMatrix1(1,:),'Parent',axes2,'LineWidth',2);