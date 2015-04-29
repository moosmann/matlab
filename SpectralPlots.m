function SpectralPlots(phase_exact,phase,y)
% Plot the line spectrum (log(abs(fft()))) of the exact and reconstructed
% phases

[dimx,dimy] = size(phase_exact);
if (nargin<3) || isempty(y),y = floor(dimy/2);end;
abscissa_length  = 60;
% number of plots in the figure
if dimx/2/abscissa_length < 4
    nmax = floor(dimx/2/abscissa_length);
else 
    nmax = 4;
end;
% Degree of polynomial fit
degree = 0;
phase_rec_1 = phase(:,:,1);
phase_rec_123 = phase(:,:,1)+phase(:,:,2)+phase(:,:,3)+phase(:,:,4);

linespec_exact   = log10(abs(fft(phase_exact(:,y))));
linespec_exact   = linespec_exact(:);
linespec_rec_0   = log10(abs(fft(phase_rec_1(:,y))));
linespec_rec_0   = linespec_rec_0(:);
linespec_rec_012 = log10(abs(fft(phase_rec_123(:,y))));
linespec_rec_012 = linespec_rec_012(:);
% Show line spectra
figure('Name','LogPlot of Line Spectra: black = exact, blue = BRO, red = BROCOR'), ...
for n = 1:nmax
    abscissa = abscissa_length*(n-1)+(1:abscissa_length);
    subplot(nmax,1,n), ... 
    plot(abscissa,linespec_exact(abscissa),'black', ...
         abscissa,linespec_rec_0(abscissa),'blue', ...
         abscissa,linespec_rec_012(abscissa),'red'),
end;
figure('Name','Errors of Line Spectra (exact - reco): blue = BRO, red = BROCOR'), ...
for n = 1:nmax
    abscissa = abscissa_length*(n-1)+(1:abscissa_length);
    subplot(nmax,1,n), ... 
    plot(abscissa,abs(linespec_exact(abscissa)-linespec_rec_0(abscissa)),'blue', ...
         abscissa,abs(linespec_exact(abscissa)-linespec_rec_012(abscissa)),'red'),
end;
% If degree=0, no fit is done nor printed.
if degree
% Fit polynomial to data set
[linespec_exact_pfit,S_exact,mu_exact] = ... 
    polyfit(1:length(linespec_exact),linespec_exact',degree);
[linespec_rec_0_pfit,S_0,mu_0] = ... 
    polyfit(1:length(linespec_rec_0),linespec_rec_0',degree);
[linespec_rec_012_pfit,S_012,mu_012] = ... 
    polyfit(1:length(linespec_rec_012),linespec_rec_012',degree);
% Show polynomial fits to line spectra
% Since the fit is shit, trigonometric functions should be fitted 
figure('Name','Polynomial Fit to Line Spectra: black = exact, blue = BRO, red = BROCOR'), ...
for n = 1:nmax
    abscissa = abscissa_length*(n-1)+(1:abscissa_length);
    subplot(nmax,1,n), ... 
    plot(abscissa,polyval(linespec_exact_pfit,abscissa),'black', ...
         abscissa,polyval(linespec_rec_0_pfit,abscissa),'blue', ...
         abscissa,polyval(linespec_rec_012_pfit,abscissa),'red'),
end;
end;


