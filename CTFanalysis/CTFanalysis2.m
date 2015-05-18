function [out, afintr, afintc, afint, int, phase] = CTFanalysis2(distanceList_m,VariableSineArg)
%Analysis of the CTF model for increasing values of the absolute phase
%shift.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    distanceList_m = 1;
    %distanceList_m = 0.25:0.25:5;
end
if nargin < 2
    VariableSineArg = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over distances
for dd = numel( distanceList_m ):-1:1
    distance_m = distanceList_m(dd);
    out(dd).distance_m = distance_m;

%% Parameters.
doFits = 1;
filestring = [ getenv('HOME'), '/data/test_pattern/lena/lena.tif' ];
blurring = [8 8 1];
energy_keV = 20;
pixelsize_m = 1e-6;
padding = 'padzeros';'padsym';'none'; 
% Read image, padding
switch lower(padding)
    case {'none',0}
        phase0 = normat( double( imread( filestring ) ) );
    case 'padzeros'
        dimX = 1024;
        phase0 = zeros(dimX,dimX);
        phase0(dimX/2+(-255:256),dimX/2+(-255:256)) = normat(double(imread(filestring)));        
    case 'padsym'
        dimX = 1024;
        phase0 = normat( double( imread( filestring ) ) );
        padSize = (dimX -size( phase0, 1) ) / 2 ;
        phase0 = padarray( phase0, [padSize padSize], 'symmetric', 'both' );        
end
dimX = size(phase0,1);
lambda = EnergyConverter(energy_keV);
SineArgPreFac = pi * lambda * distance_m / pixelsize_m^2; % normalised frequencies ranging from -1/2 to 1/2
pixPosZeroCross = @(n) round( dimX * sqrt(n* pi / SineArgPreFac ) );
% Maximum phase shifts
%rescalVec = [1, 10:10:501];
rescalVec = 1:500;
out(dd).S = rescalVec;
MaxPhaseShift = 0.01*rescalVec;
numMaxPhaseShift = numel( MaxPhaseShift );
out(dd).MaxPhaseShift = MaxPhaseShift;
% Save path
OutputPath = '/home/jmoosmann/data/QP';
switch VariableSineArg
    case 0
        coefStr = '_SinArgCoefFix';
    case 1
        coefStr = '_SinArgCoefVar';
end
OutputPath = sprintf('%s/E%ukev_z%03.0fcm_px%0.1fmu%s',OutputPath,energy_keV, distance_m * 100, pixelsize_m*1e6, coefStr );
OutputPath = regexprep(OutputPath,'\.','p');
CheckAndMakePath(OutputPath);
out(dd).OutputPath = OutputPath;


%% Blur image.
if blurring(1) > 0
    hsizex = blurring(1);
    hsizey = blurring(2);
    sigma  = blurring(3);
    phase  = imfilter(phase0,fspecial('gaussian',[hsizex hsizey],sigma),'replicate');
    phase = normat(phase);
end;

%% Print info
fprintf( '\nPARAMETERS:\n Energy: %g keV\n Distance: %g m\n Pixelsize: %g m', energy_keV, distance_m, pixelsize_m )
fprintf( '\n Padded image size: [%g %g] pixel', size(phase) )
fprintf( '\n Maximum phase shift rad: %g', max( MaxPhaseShift ) )
fprintf( '\n Number of intensity maps to compute: %g', numMaxPhaseShift )
fprintf( '\n Prefactor in argument of sine: %g', SineArgPreFac )
fprintf( '\n Number of zero crossing along axis: %g', floor( SineArgPreFac / pi /4 ) )
fprintf( '\n Pixel position of first zero crossing: %g', pixPosZeroCross(1) )
fprintf( '\n Absolute phase shift: [Min Max] = [%g %g]', MaxPhaseShift(1), MaxPhaseShift(end) )
fprintf( '\n Output path: %s', OutputPath )
%% Intensities
tic
int = zeros( dimX, dimX, numMaxPhaseShift );
fprintf('\nFORWARD PROPAGATION.')
parfor kk = 1:numMaxPhaseShift
    % Forward propagation
    int(:,:,kk) = Propagation( MaxPhaseShift(kk) * phase, [energy_keV distance_m pixelsize_m], 1, 'symmetric', 0 );
    % Fourier transform
end
fprintf('\n Time elapsed: %g s',toc)

%% FT intensity
afint = zeros( dimX, dimX, numMaxPhaseShift );
for kk = 1:numMaxPhaseShift
    % Modulus of Fourier transform of intensity contrast
    afint(:,:,kk) = (abs( fftshift( fft2( SubtractMean( int(:,:,kk) ) ) ) ) );    
end
fprintf('\nFOURIER TRANSFORM.\n Time elapsed: %g s',toc)

%% Angular integration
fprintf('\nANGULAR INTEGRATION.')
radiusStepSize = 1;
numThetaSteps = round( pi * dimX /1 );
thetaoffset = pi * 5 /180;
thetastep = (pi/2 - 2 * thetaoffset) / (numThetaSteps/4);
theta = thetaoffset : thetastep : pi/2-thetaoffset;
theta = cat( 2, theta, theta+pi/2, theta+pi, theta+3*pi/2 );
numThetaSteps = numel(theta);
radiusMax = round( dimX/2 - 1 );
radius = 1 : radiusStepSize : radiusMax;
numRadialSteps = numel(radius);
% alocate
afintc = zeros( numRadialSteps, numMaxPhaseShift, numThetaSteps);
fprintf('\n Number of angular steps: %u',numThetaSteps)
% Start looping
parfor kk = 1:numMaxPhaseShift;
    for rr = 1:numRadialSteps
        afintkk = afint( :, :, kk);
        radiusrr = radius(rr);
        for tt = 1:numThetaSteps
            thetatt = theta(tt);
            %thetatt = 2 * pi * tt / numThetaSteps; 
            afintc(rr,kk,tt) = afintkk( round(dimX/2+radiusrr*cos(thetatt)), round(dimX/2+radiusrr*sin(thetatt)));
        end
        %afintr(rr,kk) = sum(afintc) / numThetaSteps;
    end
end
% Integrate
afintr = sum( afintc, 3) / numThetaSteps;
% Normalise
afintr = afintr ./ repmat( max( afintr, [], 1), [size(afintr,1) 1] );
fprintf('\n Time elapsed: %g s',toc)

%% Fit integrated line.
if doFits    
    x = pixPosZeroCross(0.4) /radiusStepSize : pixPosZeroCross(1.75) /radiusStepSize;
    xMinSearch = pixPosZeroCross(0.5) /radiusStepSize : pixPosZeroCross(1.4) /radiusStepSize;
    out(dd).x = x;
    xn = SineArgPreFac * ( radiusStepSize * x / dimX ) .^2 / pi;
    xn = xn(:);
    out(dd).xn = xn;
    out(dd).pixelPosToSinArg = @(x) SineArgPreFac * ( radiusStepSize * x/ dimX ) .^2 / pi;
    out(dd).pixelPosToSinArgROI = @(x) SineArgPreFac * ( radiusStepSize * (x + xMinSearch(1)-1)/ dimX ) .^2 / pi;
    xnMinSearch = SineArgPreFac * ( radiusStepSize * xMinSearch / dimX ) .^2 / pi;    
    fprintf('\nFIT INTEGRATED LINE')
    fprintf('\n Support of fit: pixel position [%u .. %u], argument [%g ..%g]',x(1),x(end),xn(1),xn(end))
    cf{numMaxPhaseShift} = [];
    parfor kk = 1:numMaxPhaseShift
        y = afintr(x,kk) ;
        %y = y / max(y(:));
        cf{kk} = FitLine( xn, y(:) , VariableSineArg);
    end
    out(dd).cf = cf;
    fprintf('\n Time elapsed: %g s',toc)
     
    % Determine values and positions aof minima and maxima.     
    for kk = numMaxPhaseShift:-1:1
        [out(dd).cfMaxVal(kk), out(dd).cfMaxPos(kk)] = max(cf{kk}(xnMinSearch));
        [out(dd).cfMinVal(kk), out(dd).cfMinPos(kk)] = min(cf{kk}(xnMinSearch));
    end
end
fprintf('\n')

for kk = 1:numMaxPhaseShift
    plot( xn, afintr(x,kk) ,'.', xn, cf{kk}(xn) ),
    axis tight equal,
    pause(0.1)    
    %saveas(gcf,sprintf('%s/plot_%03u.png', OutputPath, kk),'png')
end

plot(out(dd).MaxPhaseShift,out(dd).cfMinPos,'.-')
saveas( gcf, sprintf( '%s/MinPos_vs_MaxPhaseShift.png', OutputPath), 'png' )
OutputPath = '/home/jmoosmann/data/QP/NormalisedMinPos_vs_MaxPhaseShift';
CheckAndMakePath( OutputPath )
plot( out(dd).MaxPhaseShift, normat( out(dd).cfMinPos ), '.-' )
saveas( gcf, sprintf( '%s/MinPos_vs_MaxPhaseShift_z%03.0fcm.png', OutputPath, distance_m * 100 ),'png')

end

%% Figure 1b
figure('Name','Figure 1b: data and fits')
ss = [1,100,370,500];
for nn = numel(ss):-1:1
    [~, Spos(nn)] = max(out.S == ss(nn));
    Y(:,nn) = out.cf{Spos(nn)}(out.xn);
end
h = plot(out.xn,afintr(out.x,Spos),'.',out.xn,Y,'-');
axis equal tight
h(1).Color='r';
h(2).Color='b';
h(3).Color='m';
h(4).Color='g';
h(5).Color='r';
h(6).Color='b';
h(7).Color='m';
h(8).Color='g';
%set(0, 'defaultTextInterpreter', 'tex');
xlabel('\lambda z {\bf \xi }^2')
legend(h,'data: S = 1','data: S = 100','data: S = 370','data: S = 500','fit: S = 1','fit: S = 100','fit: S = 370','fit: S = 500')
OutputPath = '/home/jmoosmann/data/QP/figures';
saveas( gcf, sprintf( '%s/Fig-1b.png', OutputPath),'png')

%% Figure 1c
figure('Name','Figure 1c: minima position')
plot(out.S,out.pixelPosToSinArgROI(out.cfMinPos),'.-');
xlabel('S')
ylabel('\lambda z |{\bf\xi}|_{1,min}^2')
saveas( gcf, sprintf( '%s/Fig-1c.png', OutputPath),'png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitresult, gof] = FitLine(x, y, VariableSineArg)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'fit' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.

%  Auto-generated by MATLAB on 07-Apr-2015 16:09:59

if nargin < 3
    VariableSineArg = 1;
end

%% Fit: 'fit'.
[xData, yData] = prepareCurveData( x, y);

% Set up fittype and options.
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
switch VariableSineArg
    case 1
        ft = fittype( 'a1 * exp( -a2 * pi*x  ) * ( abs( sin(c0*pi*x) )  + b0 + b1 * sqrt(pi*x) + b2 * pi*x + b3 * sqrt(pi*x)^3 + b4 * (pi*x)^2 + b5  * sqrt(pi*x)^5 )', 'independent', 'x', 'dependent', 'y' );
        opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 0];
        opts.StartPoint = [1 1 1 0 0 0 0 0 0.999];       
        opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf 1];
    case 0
        ft = fittype( 'a1 * exp( -a2 * pi*x  ) * ( abs( sin(   pi*x) )  + b0 + b1 * sqrt(pi*x) + b2 * pi*x + b3 * sqrt(pi*x)^3 + b4 * (pi*x)^2 + b5  * sqrt(pi*x)^5 )', 'independent', 'x', 'dependent', 'y' );
        opts.StartPoint = [1 1 1 0 0 0 0 0];        
end

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

