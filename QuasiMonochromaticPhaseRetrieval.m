ca;clear all
% directory and folder
ParPath = '/mnt/tomoraid-LSDF/users/moosmann/phd/figures/QuasiInt/';
% energy range
Emean = 20; % kev
bandwidth = 1e-1;
bwstr = regexprep(sprintf('%f',bandwidth),'\.','p');
Efwhm = bandwidth * Emean;
Esigma = Efwhm/(2*sqrt(2*log(2)));

%% Refractive index
% Use values of refractive index of water from Henke
data = importdata([ParPath 'refractive_index.dat'],' ',2);
delta = data.data(:,2);
% normalize delta to value of delta at E_mean
delta  = delta/delta(501);

%% Energy spectrum
% Weight function: Gaussian
gp = @(x,x0,sig) 1/(sig*sqrt(2*pi))*exp(-((x-x0).^2)/(2*sig^2));
% Range of energy including nsigma times the standar deviation
nsigma = 5;
N = 500;
dx = 2*nsigma*Esigma/N;
x = (Emean - nsigma*Esigma):dx:(Emean + nsigma*Esigma);
dE = x(2)-x(1);
% Gaussian weights
y = gp(x,Emean,Esigma);
% plot
%plot(x,y)
fprintf('Energy, mean: %g keV\n',Emean)
fprintf('Energy, FWHM: %g keV\n',Efwhm)
fprintf('Energy, sigma: %g keV\n',Esigma)
fprintf('Energy, [min max]: [%g %g] keV\n',x(1),x(end))
fprintf('Energy, step size: %g keV (%g keV) keV\n',dE,dx)
fprintf('Number of points: %u, Number of sigma: %u \n',N,nsigma)
fprintf('Area under Gaussian: %g\n',sum(y)*dx)
%% Quasimonochromatic intensity
MaxPhaseShift = 0.1;
distance = 1; %m
pixelsize = 1e-6; %m
pha0 = MaxPhaseShift*normat(double(imread('~/lena.tif')));
intMean = 0;
for nn = numel(x):-1:1
    int(:,:,nn) = Propagation(delta(nn)*pha0,[x(nn) distance pixelsize],2,'symmetric',0);
    intMean = intMean + int(:,:,nn)*y(nn)*dx;
end
%intMean = intMean*dx;
intEmean = int(:,:,ceil((N+1)/2));
intDiff = 100*abs(intEmean-intMean)./intEmean;

% Plot
itool([intEmean intMean]),
itool(intDiff)

%% save images and colorbar

% image and colorbar
dynRange(1) = min([intEmean(:); intMean(:)]);
dynRange(2) = max([intEmean(:); intMean(:)]);

imwrite(normat(intEmean,dynRange),[ParPath 'IntMono.png'])
imwrite(normat(intMean,dynRange), sprintf('IntPoly_%s.png',ParPath,bwstr))
savecolorbar(dynRange,[ParPath 'Int_Colorbar']);
% difference map
imwrite(normat(intDiff),[ParPath 'Diff_IntMono_IntPoly.png'])
savecolorbar([min(intDiff(:)) max(intDiff(:))],[ParPath 'Diff_Colorbar']);


%% phase retrieval
pf = PhaseFilter('ctf',size(intMean),[Emean distance pixelsize],2.5,0.1);
phaMean = ifft2(pf.*fft2(intMean));
phaEmean = ifft2(pf.*fft2(intEmean));
%itool(phaMean),itool(phaEmean)

if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BM spectrum
    data = importdata([ParPath 'Spectrum_BM_APS_E8to130kev_steps500.dat'],' ',2);
    x2 = data.data(:,1)/1000;% E in keV
    dx2 = x(2)-x(1);
    y2 = data.data(:,2);
    % normalize spectrum
    y2 = y2/sum(y2);
    % intensity
    intMean2 = 0;
    for nn = numel(x2):-1:1
        int2(:,:,nn) = Propagation(pha0,[x2(nn) distance pixelsize],2,'symmetric',0);
        intMean2 = intMean2 + int2(:,:,nn)*y2(nn);
    end
end