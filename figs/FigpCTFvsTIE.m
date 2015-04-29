function FigpCTFvsTIE(EnergyDistancePixelsize,dimx)

%% Default arguments.
if nargin < 1
    EnergyDistancePixelsize = [22 500e-3 1.5e-6];
end
if nargin < 2
    dimx = 20*2048;
end

ca
%% Parameters
xcut      = 1;
alphaTIE  = 2.5;
alphaCTF  = alphaTIE;
BinFiltThres = 0.3;
% Wave length.
lambda    = 6.62606896e-34*299792458/(EnergyDistancePixelsize(1)*1.60217733e-16);
% Prefactor needed for TIE and CTF retrieval.
prefactor = EnergyDistancePixelsize(3)^2/(2*pi*lambda*EnergyDistancePixelsize(2));

%% Filters.
%x       = (-1/2:1/dimx:1/2-1/dimx);
x       = 1/sqrt(prefactor)*(1/dimx:1/dimx:1/2-1/dimx);
x2      = x.^2/2;
xi2     = (2*x2 + 10^-alphaTIE);
invxi2  = 1./xi2;
%invsine =1./(2*sign(sin(1/prefactor*x.^2/2)).*(abs(sin(1/prefactor*x.^2/2))) + 10^-alphaCTF);
sine    = ((2*sign(sin(x2)).*(abs(sin(x2))) + 10^-alphaCTF));
asine   = abs((2*(sign(sin(x2))).*(abs(sin(x2))) + 10^-alphaCTF));
BinFilt = ones(1,length(x));
BinFilt( (sin(x2).^2<=BinFiltThres) & ((x)>pi/2) ) = 0;
invsine = 1./sine;

%% Plot.
figure1 = figure('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98404194812 29.67743169791],...
    'Name','FT intensity modulation according TIE, CTF, pCTF');
axes1 = axes('Parent',figure1,'FontSize',12);
xx = 1:floor(length(x)/xcut);
plot1 = plot(x(xx),xi2(xx),'black',x(xx),1.01*asine(xx),'blue',x(xx),BinFilt(xx).*asine(xx),'red','Parent',axes1);
set(plot1,'LineWidth',2)
xlabel('sqrt(\pi \lambda z {\bf k}^2)','FontSize',20);
%ylabel('y','FontSize',18);
ylim(axes1,[0 4]);
legend('linearized TIE','CTF','projected CTF','Location','NorthEast')
%set(gcf,'defaulttextinterpreter','none')
saveas(gcf,'/home/moosmann/seminartalk/FourierModulation_quadr.eps','epsc2')
set(gcf,'PaperPositionMode','Auto') 
saveas(gcf,'/home/moosmann/seminartalk/FourierModuation_nonquadr.eps','epsc2')

figure2 = figure('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98404194812 29.67743169791],...
    'Name','Fourier space filters according TIE, CTF, pCTF applied to FT intensity');
axes2 = axes('Parent',figure2,'FontSize',12);
xx    = 1:floor(length(x)/xcut);
plot2 = plot(x(xx),invxi2(xx),'black',x(xx),invsine(xx),'blue',x(xx),BinFilt.*invsine(xx),'red','Parent',axes2);
set(plot2,'LineWidth',2,'LineStyle','-')
set(plot2(2),'LineWidth',1.5)
set(plot2(2),'LineStyle','none','Marker','.','MarkerSize',5)
set(plot2(3),'LineStyle','none','Marker','.','MarkerSize',5)
xlabel('sqrt(\pi \lambda z {\bf k}^2)','FontSize',20);
%ylabel('y','FontSize',18);
ylim(axes2,[-3 3]);
legend('linearized TIE','CTF','projected CTF','Location','SouthWest')
%set(gcf,'defaulttextinterpreter','none')
saveas(gcf,'/home/moosmann/seminartalk/FourierFilters_quadr.eps','epsc2')
set(gcf,'PaperPositionMode','Auto') 
saveas(gcf,'/home/moosmann/seminartalk/FourierFilters_nonquadr.eps','epsc2')


% set(gca,'XTick',0:pi/4:pi)
% set(gca,'XTickLabel',{'0','pi/4','pi','3pi/4','pi'})
% for x=[pi/12 pi/6 pi/4 pi/2]
%        text(x,(x-sin(x))/x,['(' num2str(x,'%.2g') ',' num2str((x-sin(x))/x,'%3.2g') ...
%                       ',' num2str((x-sin(x))/sin(x),'%3.2g') ')']);
% end;
 