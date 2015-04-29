
%figure,plot(x,1./(1+exp(-2*k*x)),x,sinh(2*k*x/2)./cosh(2*k*x/2).^3);
clear all
% Poisson noise level
NoiseLevel = 1000;
noiseim = @(im,NoiseLevel) 1/NoiseLevel*double(imnoise(uint16(NoiseLevel*im),'poisson'));
% Phantom
switch 'edge'
    case 'edge'
        % Size of images
        dim1 = 512;
        dim2 = 2048;
        xhalf = round(dim2/2);
        yhalf = round(dim1/2);
        [x,y] = meshgrid((1:dim2)-xhalf,1:dim1);
        %% phantom
        % Model an edge-like transition by a Fermi-Dirac-Distribution: H_k(x) =
        % 1./(1+exp(-2*k*x), where k defines the sharpness of the edge. Set k such
        % that a distance dx away from the edge the value of the edge has dropped
        % by a factor of dH
        k = @(dx,dH) 1/(2*dx)*log(1/dH-1);
        fdedge = @(phimax,dx,dH) phimax*1./(1+exp(-2*k(dx,dH)*x));
        phan1 = fdedge(1,4,1/25);
        phan2 = fdedge(1,16,1/16);
    case 'lena'
        % Size of images
        dim1 = 512;
        dim2 = 512;
        xhalf = round(dim2/2);
        yhalf = round(dim1/2);
        [x,y] = meshgrid((1:dim2)-xhalf,1:dim1);
        phan1 = imread('~/lena.tif');
        phan1 = phan1(yhalf,:);
        phan1 = normat(double(repmat(phan1,[dim1 1])));
end
if 0
    xplot = xhalf+(-199:200);
    figplot('Phantom line cuts',[phan1(yhalf,xplot)' phan2(yhalf,xplot)'])
end
%% intensity
z0 = 0.1;
zMax = 3;
m = 100;
z = z0 + (zMax - z0)*(0:m-1)/(m-1);
int1 = zeros(dim1,dim2,m);
nint1 = zeros(dim1,dim2,m);
noise1 = zeros(dim1,dim2,m);
fnoise1 = zeros(dim1,dim2,m);
fint1 = zeros(dim1,dim2,m);
fnint1 = zeros(dim1,dim2,m);
parfor nn = 1:m
    int1(:,:,nn) = Propagation(phan1,[20 z(nn) 1e-6],2,'symmetric',0);
    nint1(:,:,nn) = noiseim(int1(:,:,nn),NoiseLevel);
    noise1(:,:,nn) = nint1(:,:,nn) - int1(:,:,nn);
    fint1(:,:,nn) = fftshift(fft(int1(:,:,nn)-1,[],2));
    fnint1(:,:,nn) = fftshift(fft(nint1(:,:,nn)-1,[],2));
    fnoise1(:,:,nn) = fftshift(fft(noise1(:,:,nn),[],2));
end
if 0
    int2 = zeros(dim1,dim2,m);
    nint2 = zeros(dim1,dim2,m);
    fint2 = zeros(dim1,dim2,m);
    fnint2 = zeros(dim1,dim2,m);
    noise2 = zeros(dim1,dim2,m);
    fnoise2 = zeros(dim1,dim2,m);
    parfor nn = 1:m
        int2(:,:,nn) = Propagation(phan2,[20 z(nn) 1e-6],2,'symmetric',0);
        nint2(:,:,nn) = noiseim(int2(:,:,nn),NoiseLevel);
        noise2(:,:,nn) = nint2(:,:,nn) - int2(:,:,nn);
        fint2(:,:,nn) = fftshift(fft(int2(:,:,nn)-1,[],2));
        fnint2(:,:,nn) = fftshift(fft(nint2(:,:,nn)-1,[],2));
        fnoise2(:,:,nn) = fftshift(fft(noise2(:,:,nn),[],2));
    end;
end
%% 2D Plots
%Phantom 1
mroi = [1 10 50 100];
xplot = xhalf + (-99:100);
X = x(yhalf,xplot);
Y = squeeze(int1(yhalf,xplot,mroi));
zroi = z(mroi);
figure('Name','Phantom 1. Intensity line cuts')
plot(X,Y)
fdistPlot = max(Y) - min(Y);
fdistTheo = sqrt(EnergyConverter(20)*zroi)*10^6;

% Phantom 1+2
if 0
    % Line cuts through maps of contrast evolution
    %nimplay(int1(:,:,:));
    intMax = max(max(int1(yhalf,xplot,:)));
    intMin = min(min(int1(yhalf,xplot,:)));
    phan1Cut = squeeze(phan1(yhalf,xplot)');
    intCut1 = 1*(squeeze(int1(yhalf,xplot,:))-1);
    intCut2 = 1*(squeeze(int2(yhalf,xplot,:))-1);
    plotIntCut = @(nn) figplot(sprintf('Intensity line cuts. z:%gm',z(nn)),[intCut1(:,nn) intCut2(:,nn)]);
    plotIntCut(1)
    plotIntCut(2)
    plotIntCut(floor(m/2))
end

%% Phantom 1: Surface plots
dxplot = floor(0.3*dim2/2);
xplot = xhalf+(-dxplot:dxplot);
[xroi,yroi] = meshgrid(z,squeeze(x(1,xplot)));
xplotf = 1:2:dim2;
[xroif,yroif] = meshgrid(z,squeeze(x(1,xplotf)));
if 1
    % NOISELESS INTENSITY
    % Fringe evolution
    figure('Name','Phantom 1. Noiseless. Surface plot of fringe evolution')
    mesh( xroi, yroi, squeeze( int1(yhalf,xplot,:) ) )
    % Evolution of spectrum
    figure('Name','Phantom 1. Noiseless. Surface plot of spectrum (abs(fft))')
    mesh(xroif,yroif,log( 1 + abs( ( squeeze( fint1(yhalf,xplotf,:) ) ) ) ))
    % POISSON NOISED INTENSITY
    % Fringe evolution
    figure('Name','Phantom 1. Noised. Surface plot of fringe evolution')
    mesh( xroi, yroi, squeeze( nint1(yhalf,xplot,:) ) )
    % Evolution of spectrum
    figure('Name','Phantom 1. Noised. Surface plot of spectrum (abs(fft))')
    mesh(xroif,yroif,log( 1 + abs( ( squeeze( fnint1(yhalf,xplotf,:) ) ) ) ))
    % POISSON NOISE
    % Noise evolution
    figure('Name','Phantom 1. Surface plot of noise evolution')
    mesh( xroi, yroi, squeeze( noise1(yhalf,xplot,:) ) )
    % Evolution of noise spectrum
    figure('Name','Phantom 1. Surface plot of noise spectrum (abs(fft))')
    mesh(xroif,yroif,log( 1 + abs( ( squeeze( fnoise1(yhalf,xplotf,:) ) ) ) ))
end

%% Phantom 2: Surface plots
if 0
    % NOISELESS INTENSITY
    % Fringe evolution
    figure('Name','Phantom 2. Noiseless. Surface plot of fringe evolution')
    mesh( xroi, yroi, squeeze( int2(yhalf,xplot,:) ) )
    % Evolution of spectrum
    figure('Name','Phantom 2. Noiseless. Surface plot of spectrum (abs(fft))')
    mesh(xroi,yroi,log( 1 + abs( ( squeeze( fint2(yhalf,xplot,:) ) ) ) ))
    % POISSON NOISED INTENSITY
    % Fringe evolution
    figure('Name','Phantom 2. Noised. Surface plot of fringe evolution')
    mesh( xroi, yroi, squeeze( nint2(yhalf,xplot,:) ) )
    % Evolution of spectrum
    figure('Name','Phantom 2. Noised. Surface plot of spectrum (abs(fft))')
    mesh(xroi,yroi,log( 1 + abs( ( squeeze( fnint2(yhalf,xplot,:) ) ) ) ))
    % POISSON NOISE
    % Noise evolution
    figure('Name','Phantom 1. Surface plot of noise evolution')
    mesh( xroi, yroi, squeeze( noise2(yhalf,xplot,:) ) )
    % Evolution of noise spectrum
    figure('Name','Phantom 1. Surface plot of noise spectrum (abs(fft))')
    mesh(xroi,yroi,log( 1 + abs( ( squeeze( fnoise2(yhalf,xplot,:) ) ) ) ))
end

%% 2D plots: Spectrum and noise
% PHANTOM 1
if 0
    % noiseless
    implot = @(m) (log(1+abs(fftshift(fft(SubtractMean((squeeze(int1(yhalf,501:1500,m)))))))));
    figure('Name','Spectrum (abs(fft)) of noiseless intensity'),
    plot1 = plot( [implot(1); implot(round(m/2)); implot(m)]');
    % set(plot1(1),'Color',[1 0 0]);
    % set(plot1(2),'Color',[0 1 0]);
    % set(plot1(3),'Color',[0 0 1]);
    % Poisson noised
    nimplot = @(m) (log(1+abs(fftshift(fft(SubtractMean(noiseim(squeeze(int1(yhalf,501:1500,m)),NoiseLevel)))))));
    figure('Name',sprintf('Spectrum (abs(fft)) of Poisson-noised intensity. Average counts: %u',NoiseLevel))
    plot2 = plot( [nimplot(1); nimplot(round(m/2)); nimplot(m)]');
    % set(plot2(1),'Color',[1 0 0]);
    % set(plot2(2),'Color',[0 1 0]);
    % set(plot2(3),'Color',[0 0 1]);
end
% PHANTOM 2
if 0
    % noiseless
    implot = @(m) (log(1+abs(fftshift(fft(SubtractMean((squeeze(int2(yhalf,501:1500,m)))))))));
    figure('Name','Spectrum (abs(fft)) of noiseless intensity'),
    plot3 = plot( [implot(1); implot(round(m/2)); implot(m)]');
    % Poisson noised
    nimplot = @(m) (log(1+abs(fftshift(fft(SubtractMean(noiseim(squeeze(int2(yhalf,501:1500,m)),NoiseLevel)))))));
    figure('Name',sprintf('Spectrum (abs(fft)) of Poisson-noised intensity. Average counts: %u',NoiseLevel))
    plot4 = plot( [nimplot(1); nimplot(round(m/2)); nimplot(m)]');
end