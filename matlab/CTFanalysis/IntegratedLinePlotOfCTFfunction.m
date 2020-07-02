% General model:
%      f(x) = a*exp(-b2*x-b1*sqrt(x))*abs(sin(w*x+p0))*(c0+c1*sqrt(x)+c2*x+c3*sqrt(x)
%                     ^3+c4*x^2+c5*sqrt(x)^5)
% Coefficients (with 95% confidence bounds):
%        a =       1.976  (-2.034e+05, 2.034e+05)
%        b1 =    3.66e-14  (fixed at bound)
%        b2 =   3.441e-14  (fixed at bound)
%        c0 =     -0.9999  (-1.029e+05, 1.029e+05)
%        c1 =       15.61  (-1.607e+06, 1.607e+06)
%        c2 =      -31.66  (-3.26e+06, 3.26e+06)
%        c3 =        27.3  (-2.81e+06, 2.81e+06)
%        c4 =      -10.61  (-1.092e+06, 1.092e+06)
%        c5 =       1.528  (-1.573e+05, 1.573e+05)
%        p0 =    -0.07791  (-0.5569, 0.4011)
%        w =       1.025  (0.867, 1.184)
% 
% Goodness of fit:
%   SSE: 140.5
%   R-square: 0.2283
%   Adjusted R-square: 0.1863
%   RMSE: 0.9777
%% Prior to evaluation the script 'CTFanlayisis.m' needs to be run or the
%% data 'CTFanalysis.mat' needs to be imported: load('CTFanalysis.mat')
clear X1 X2 YarMatrix1 plot plot1 axes1 figure1
X1=x;
SineArgPreFac = pi*EnergyConverter(energy)*distance/(paddim*pixelsize)^2;
figure1 = figure('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98404194812 29.67743169791],...
    'Name','Angular integrated line plot of FT of intensity contrast VS frequency');
% Create axes
X2 = SineArgPreFac*X1.^2;
% axes1 = axes('Parent',figure1, ...
%      'YTickLabel',{'0','200','400','600','800','1000','y'},...
%      'XTick',[25 50 75 100 125 150 175], ...
%      'FontSize',18);
axes1 = axes('Parent',figure1,'FontSize',18);
% xlim(axes1,[25 180]);
% ylim(axes1,[0 1200]);
box(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
%YarMatrix1(6,:) = cfar{200}(X1);
%YarMatrix1(5,:) = cfar{100}(X1);
%YarMatrix1(4,:) = cfar{1}(X1);
YarMatrix1(3,:) = yar(X1,200);
YarMatrix1(2,:) = yar(X1,100);
YarMatrix1(1,:) = yar(X1,1);
% Some options
set(gcf,'defaulttextinterpreter','none')
set(gcf,'position',[0 0, 1600 900]) 
% Create multiple lines using matrix input to plot
plot1 = plot(X2,YarMatrix1,'Parent',axes1,'LineWidth',2);
set(plot1(1),'MarkerSize',15,'Marker','.','LineStyle','none','Color',[0 0 0],'LineWidth',0.5);
set(plot1(2),'MarkerSize',8,'Marker','x','LineStyle','none','Color',[0 0 1],'LineWidth',0.5);
set(plot1(3),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 0 0],'LineWidth',0.5);
%set(plot1(4),'LineStyle','--','Color',[1 0 0]);
%set(plot1(5),'LineStyle','-.','Color',[0 1 0]);
%set(plot1(6),'LineStyle','-.','Color',[0 0 1]);
% Create xlabel
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
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
% Save figure
saveas(gcf,'/home/moosmann/oe2011HMB/Fig-2new_quadr.eps','epsc2')
set(gcf,'PaperPositionMode','Auto') 
saveas(gcf,'/home/moosmann/oe2011HMB/Fig-2new.eps','epsc2')
saveas(gcf,'/home/moosmann/oe2011HMB/IntegratedLinePlotOfAbsCTF_(Fig-2new).eps','epsc2')
