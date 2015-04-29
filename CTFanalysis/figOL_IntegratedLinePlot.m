%Prior to evaluation the script 'CTFanlayisis.m' needs to be run or the
%data 'CTFanalysis.mat' needs to be imported: load('CTFanalysis.mat')

ca
clear X1 YMatrix1 plot1 axes1 figure1 X2 YarMatrix1 plot2 axes2 figure2 plot3 axes3 figure3
X1=x;
SineArgPreFac = pi*EnergyConverter(energy)*distance/(paddim*pixelsize)^2;
figure1 = figure('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98404194812 29.67743169791],...
    'Name','Angular integrated line plot of FT of intensity contrast VS frequency');

% Create axes
X2 = SineArgPreFac*X1.^2;
axes1 = axes('Parent',figure1, ...
     'YTickLabel',{'0','200','400','600','800','1000','y'},...
     'XTick',[25 50 75 100 125 150 175], ...
     'FontSize',18);
xlim(axes1,[25 180]);
ylim(axes1,[0 1200]);
box(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
YMatrix1(5,:) = 100*abs(sin(SineArgPreFac*X1.^2));
YMatrix1(1,:) = y(X1,450);
YMatrix1(2,:) = y(X1,200);
YMatrix1(3,:) = cf{450}(X1);
YMatrix1(4,:) = cf{200}(X1);
% Some options
set(gcf,'defaulttextinterpreter','none')
set(gcf,'position',[0 0, 1600 900]) 
% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',2);
set(plot1(1),'MarkerSize',8,'Marker','x','LineStyle','none','Color',[1 0 0],...
    'LineWidth',0.5);
set(plot1(2),'Marker','.','LineStyle','none','Color',[0 0 1],...
    'LineWidth',0.5);
set(plot1(3),'LineStyle','--','Color',[1 0 0]);
set(plot1(4),'LineStyle','-.','Color',[0 0 1]);
set(plot1(5),'Color',[0 0 0]);
% Create xlabel
xlabel('x','FontSize',18);
% Create arrow
annotation(figure1,'arrow',[0.765476190476191 0.803571428571429],...
    [0.457260556127703 0.247167868177137]);
% Create textarrow
annotation(figure1,'textarrow',[0.670833333333333 0.718452380952381],...
    [0.699917431192661 0.116374871266735],'TextEdgeColor','none','FontSize',18,...
    'String',{'x1'});
% Create textarrow
annotation(figure1,'textarrow',[0.316666666666667 0.282142857142857],...
    [0.80126570545829 0.630526315789474],'TextEdgeColor','none','FontSize',18,...
    'String',{'x3'});
% Create arrow
annotation(figure1,'arrow',[0.317261904761905 0.33452380952381],...
    [0.799 0.382105263157895]);
% Create textarrow
annotation(figure1,'textarrow',[0.764285714285714 0.72202380952381],...
    [0.460350154479918 0.151578947368421],'TextEdgeColor','none','FontSize',18,...
    'String',{'x2'});
% Save figure
saveas(gcf,'/home/moosmann/oe2011HMB/Fig-1b.eps','epsc2')
set(gcf,'PaperPositionMode','Auto') 
saveas(gcf,'/home/moosmann/oe2011HMB/Fig-1b_nonquadr.eps','epsc2')
saveas(gcf,'/home/moosmann/oe2011HMB/IntegratedLinePlotOfAbsFTg_(Fig-1b).eps','epsc2')

%% Plot ratio
figure2 = figure('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'Name','Angular integrated line plot of ratio of FT of intensity contrast and FT of phase VS frequency');
%'PaperSize',[20.98404194812 29.67743169791],...
axes2 = axes('Parent',figure2, ...
     'XTick',[25 50 75 100 125 150 175], ...
     'FontSize',18);
 %'YTickLabel',{'0','2','4','y'},...
xlim(axes2,[25 180]);
ylim(axes2,[0 6]);
box(axes2,'on');
hold(axes2,'all');
% Create multiple lines using matrix input to plot
YarMatrix1(5,:) = 2*abs(sin(X2));
YarMatrix1(1,:) = yar(X1,450);
YarMatrix1(2,:) = yar(X1,200);
YarMatrix1(3,:) = cfar{450}(X1);
YarMatrix1(4,:) = cfar{200}(X1);
% Some options
set(gcf,'defaulttextinterpreter','none')
set(gcf,'position',[0 0, 1600 900]) 
% Create multiple lines using matrix input to plot
plot2 = plot(X1,YarMatrix1,'Parent',axes2,'LineWidth',2);
set(plot2(1),'MarkerSize',8,'Marker','x','LineStyle','none','Color',[1 0 0],...
    'LineWidth',0.5);
set(plot2(2),'Marker','.','LineStyle','none','Color',[0 0 1],...
    'LineWidth',0.5);
set(plot2(3),'LineStyle','--','Color',[1 0 0]);
set(plot2(4),'LineStyle','-.','Color',[0 0 1]);
set(plot2(5),'Color',[0 0 0]);
% Create xlabel
xlabel('x','FontSize',18);
% % Create arrow
% annotation(figure1,'arrow',[0.765476190476191 0.803571428571429],...
%     [0.457260556127703 0.247167868177137]);
% % Create textarrow
% annotation(figure1,'textarrow',[0.670833333333333 0.718452380952381],...
%     [0.699917431192661 0.116374871266735],'TextEdgeColor','none','FontSize',18,...
%     'String',{'x1'});
% % Create textarrow
% annotation(figure1,'textarrow',[0.316666666666667 0.282142857142857],...
%     [0.80126570545829 0.630526315789474],'TextEdgeColor','none','FontSize',18,...
%     'String',{'x3'});
% % Create arrow
% annotation(figure1,'arrow',[0.317261904761905 0.33452380952381],...
%     [0.799 0.382105263157895]);
% % Create textarrow
% annotation(figure1,'textarrow',[0.764285714285714 0.72202380952381],...
%     [0.460350154479918 0.151578947368421],'TextEdgeColor','none','FontSize',18,...
%     'String',{'x2'});
% % Save figure
% saveas(gcf,'/home/moosmann/oe2011HMB/Fig-1c.eps','epsc2')
% set(gcf,'PaperPositionMode','Auto') 
% saveas(gcf,'/home/moosmann/oe2011HMB/Fig-1c_nonquadr.eps','epsc2')
% saveas(gcf,'/home/moosmann/oe2011HMB/IntegratedLinePlotOfAbsFTg_(Fig-1b).eps','epsc2')

%% Plot ratio
figure3 = figure('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'Name','Anular integrated line plot of ratio of FT of intensity contrast and FT of phase VS frequency^2');
%'PaperSize',[20.98404194812 29.67743169791],...
axes3 = axes('Parent',figure3, ...
     'FontSize',18);
 %'XTick',[25 50 75 100 125 150 175], ...
 %'YTickLabel',{'0','2','4','y'},...
xlim(axes3,[0 5]);
ylim(axes3,[0 6]);
box(axes3,'on');
hold(axes3,'all');
% Create multiple lines using matrix input to plot
YarMatrix1(5,:) = 2*abs(sin(X2));
YarMatrix1(1,:) = yar(X1,450);
YarMatrix1(2,:) = yar(X1,200);
YarMatrix1(3,:) = cfar{450}(X1);
YarMatrix1(4,:) = cfar{200}(X1);
% Some options
set(gcf,'defaulttextinterpreter','none')
set(gcf,'position',[0 0, 1600 900]) 
% Create multiple lines using matrix input to plot
plot3 = plot(X2,YarMatrix1,'Parent',axes3,'LineWidth',2);
set(plot3(1),'MarkerSize',8,'Marker','x','LineStyle','none','Color',[1 0 0],...
    'LineWidth',0.5);
set(plot3(2),'Marker','.','LineStyle','none','Color',[0 0 1],...
    'LineWidth',0.5);
set(plot3(3),'LineStyle','--','Color',[1 0 0]);
set(plot3(4),'LineStyle','-.','Color',[0 0 1]);
set(plot3(5),'Color',[0 0 0]);
% Create xlabel
xlabel('x','FontSize',18);