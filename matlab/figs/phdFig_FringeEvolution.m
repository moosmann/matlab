
clear all

% File and folder names
FigName = 'FringeEvolution';
SavePath = sprintf('/mnt/tomoraid-LSDF/users/moosmann/phd/figures/%s/',FigName);
MakePath(SavePath)
% parameter
fontSize = 20;
% Phantom
% Size of images
dim1 = 64;
dim2 = 2*2048;
xhalf = round(dim2/2);
yhalf = round(dim1/2);
[x,y] = meshgrid((1:dim2)-xhalf,1:dim1);

%% phantom
% Model an edge-like transition by a Fermi-Dirac-Distribution: H_k(x) =
% 1./(1+exp(-k*x), where k defines the sharpness of the edge. Set k such
% that a distance dx away from the edge the value of the edge has dropped
% by a factor of dH
k = @(dx,dH) 1/(dx)*log(1/dH-1);
fdedge = @(phimax,dx,dH) phimax*1./(1+exp(-k(dx,dH)*x));
phan1 = fdedge(1,4,1/25);
% Poisson noise level
NoiseLevel = 1000;
noiseim = @(im,NoiseLevel) 1/NoiseLevel*double(imnoise(uint16(NoiseLevel*im),'poisson'));

%% intensity
z = [0.1 0.5 1.5 3];
m = length(z);
int1 = zeros(dim1,dim2,m);
% nint1 = zeros(dim1,dim2,m);
% noise1 = zeros(dim1,dim2,m);
% fnoise1 = zeros(dim1,dim2,m);
% fint1 = zeros(dim1,dim2,m);
% fnint1 = zeros(dim1,dim2,m);
parfor nn = 1:m
     int1(:,:,nn) = Propagation(phan1,[20 z(nn) 1e-6],2,'symmetric',0);
%     nint1(:,:,nn) = noiseim(int1(:,:,nn),NoiseLevel);
%     noise1(:,:,nn) = nint1(:,:,nn) - int1(:,:,nn);
%     fint1(:,:,nn) = fftshift(fft(int1(:,:,nn)-1,[],2));
%     fnint1(:,:,nn) = fftshift(fft(nint1(:,:,nn)-1,[],2));
%     fnoise1(:,:,nn) = fftshift(fft(noise1(:,:,nn),[],2));
end

%% 2D Plots
% Create data
xplot = xhalf + (-49:49);
X = x(yhalf,xplot);
Y1 = squeeze(int1(yhalf,xplot,:));
Y2 = phan1(yhalf,xplot)';%+max(phan1(:))/2;
% Create figure
figure1 = figure;
% ('XVisual','0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%     'PaperType','A0',...
%     'Name','PhD thesis. Figure: FringeEvolution of edge. Intensity line cuts');
% %'PaperSize',[20.98404194812 29.67743169791],...
% Create axes
axes1 = axes('Parent',figure1,'FontSize',fontSize,'XMinorTick','off','XTick',[ -40 -30 -20 -10 0 10 20 30 40 50]);
box(axes1,'on');
hold(axes1,'all');
% Create plot
%plot1 = plot(X,Y,'Parent',axes1);
[plot1 h1 h2] = plotyy(X,Y2,X,Y1,'Parent',axes1);
% Plot line width
set(h1,'LineWidth',1.5);
set(h2,'LineWidth',1.5);
% Fringe distance
[Y1max Y1maxpos] = max(Y1);
[Y1min Y1minpos] = min(Y1);
Yheight = Y1max - Y1min;
Ywidth = Y1maxpos - Y1minpos;
YwidthTheo = sqrt(EnergyConverter(20)*z)*1e6;
%figure(2),plotyy(sqrt(z),[Ywidth; YwidthTheo]',sqrt(z),Yheight)
% Figure properties
xlabel('xlabel','FontSize',fontSize)
ylabel(plot1(1),'ylabel','FontSize',fontSize)
ylabel(plot1(2),'ylabel2','FontSize',fontSize)
ylim(plot1(2),[0.99*min(Y1(:)) 1.01*max(Y1(:))]);
%legend2 = legend(arrayfun(@(val) sprintf('z = %g m',val),z,'UniformOutput',0));
% Create legend
legend(axes1,'show');
legend('location','NorthWest');
legend('boxoff');
% Create colorbar
set(axes1, 'Box', 'off');
% Turn off second x-axis labels and ticks
set(plot1(2),'xticklab',[],'xtick',[])
set(plot1,'FontSize',fontSize)
%colorbar('peer',axes1);
saveas(figure1,sprintf('%s%s',SavePath,'Fig1.eps'),'epsc2')
saveas(gcf,sprintf('%s%s',SavePath,'Fig1_gcf.eps'),'epsc2')
%set(gcf,'PaperPositionMode','Auto') 
%saveas(gcf,sprintf('%s%s%s',SavePath,FigName,'_nonquadr.eps'),'epsc2')

