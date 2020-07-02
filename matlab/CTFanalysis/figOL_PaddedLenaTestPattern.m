%Prior to execution the scrip 'FitCritExpoScript.m' needs to be run or the
%data 'CTFanlaysis2.mat' needs to be loaded.

ca
% Create figure
figure1 = figure('XVisual',...
    '0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)', ...
    'Name','Padded Lena test pattern');
colormap('gray');
%set(gcf,'position',[0 0, 900 900])
set(gcf,'defaulttextinterpreter','none')

% Create axes
axes1 = axes('Parent',figure1,'YTick',[0 200 400 600 800 1000],...
    'YAxisLocation','right',...
    'YDir','reverse',...
    'XTick',[0 200 400 600 800 1000],...
    'TickDir','out',...
    'Layer','top',...
    'FontSize',18,...
    'DataAspectRatio',[1 1 1],...
    'CLim',[0 0.01]);

xlim(axes1,[1 1024]);
ylim(axes1,[1 1024]);
box(axes1,'on');
hold(axes1,'all');

% Create image
image(phase,'Parent',axes1,'CDataMapping','scaled');

saveas(gcf,'/home/moosmann/ol/Fig-1a.eps','epsc2')
set(gcf,'PaperPositionMode','Auto') 
saveas(gcf,'/home/moosmann/ol/Fig-1a_nonquadr.eps','epsc2')
saveas(gcf,'/home/moosmann/ol/LenaTestPatternPadded_(Fig-1a).eps','epsc2')