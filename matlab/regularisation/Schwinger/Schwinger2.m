function phi = Schwinger2(int, Method, EnergyDistancePixelsize, BinaryFilterThreshold, loopOverRegPar)
% Another try to solve the regularization problem using the Schwinger
% inegral to circumvent the singularity of 1/x^2 at x = 0.
%
% Written by Julian Moosmann, first version: 2013-12-19, last version: 2013-12-19
%
% Schwinger2()

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    Method = 'tie';    
end
if nargin < 3
    EnergyDistancePixelsize = [30e3 0.7 1e-6];
end
if nargin < 4
    BinaryFilterThreshold = 0.05;
end
if nargin < 1
    filename = '~/barbara.png';
    if exist( filename, 'var' )
        int = intread(filename);
    else
        int = phantom( 'Modified Shepp-Logan', 512 );
    end        
    phi0 = 0.1*normat(single( int ));
    int = Propagation2(phi0,0,EnergyDistancePixelsize,0,1);
    phi0 = SubtractMean(phi0);
end
if nargin < 5
    loopOverRegPar = 1;
end
if nargin < 6
    show_updates = 100;
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FT of intensity
int = fft2(int-1);


% reference retrieval
if loopOverRegPar(1)
    h = @(x) x(:);
    regPar = 0:0.25:8;
    for mm = numel(regPar):-1:1
        phiReg(:,:,mm) = ifft2( PhaseFilter(Method,size(int),EnergyDistancePixelsize,regPar(mm),0.1,'double') .* int );
        er(mm) = sum( abs( h( phiReg(:,:,mm) )  - phi0(:) ) );
        er2(mm) = sum( abs( h( normat(phiReg(:,:,mm) ) ) - h( normat(phi0) ) ) );
    end
    [erMinVal, erMinPos] = min(er(:));
    fprintf('\n Mininum of error map of standard regularized retrieval at index %u, regPar %g, (and error value %g)',erMinPos,regPar(erMinPos),erMinVal)
    [er2MinVal, er2MinPos] = min(er2(:));
    fprintf('\n Mininum of error map of renormed standard regularized retrieval at index %u, regPar %g, (and error value %g)\n',er2MinPos,regPar(er2MinPos),er2MinVal)
    phiRegOpt = phiReg(:,:,er2MinPos);
else
    phiRegOpt = ifft2( PhaseFilter(Method,size(int),EnergyDistancePixelsize,4.4,BinaryFilterThreshold,'double') .* int );    
end

% create filters and masks for Schwinger integration
[signMask, pfs, imMask] = SchwingerFilter(Method,size(int),EnergyDistancePixelsize,BinaryFilterThreshold);

% mask for QP
int = imMask .* int;

% determine tMax where the exponential weight drops below numerical
% precision, the smallest values of  pfs=sqrt(abs(sin(pi lambda xi^2 /
% dx^2))) is close the value of pfs(2,1)=pfs(1,2), pfs(1,1) = 0. There
% are/could be smaller values close to higher order of zeros, but the order
% of magnitude should be the same
tMax = -1 / 2 * log( eps('single'))/pfs(2,1); % Reduce by another 1/2
%dt = 0.2;0.25;%tMax/1000;
num_iter_max = 10000;
dt = tMax / num_iter_max;
t = 0:dt:tMax;
num_t = numel( t );
% integrate using a loop
phi = 0;
phiuo = 0;
fprintf( '\n Schwinger integration parameter: [dt, tMax, num_t] = [%g, %g %g]', dt, tMax, num_t )
%phi = zeros([size(int) numel(t)]);
nn_break = 0;

figure( 'Name', 'Integration weight band' )
pfs_max = max2( pfs );
pfs_med = median( pfs(:));
pfs_mean = mean( pfs(:));
pfs_min = pfs(1,2);
wmin = exp( -t .* pfs_max)';
wmed = exp( -t .* pfs_med)';
wmean = exp( -t .* pfs_mean)';
wmax = exp( -t .* pfs_min)';
plot( [wmin, wmed, wmean, wmax ] )
title( sprintf( 'sampling points: %u', numel( t ) ))
legend( {sprintf( 'low: %f', pfs_max ), sprintf('med: %g', pfs_med ), sprintf('mean: %g', pfs_mean ), sprintf('high: %g', pfs_min )} )

figure( 'Name', 'Schwinger iteration' )
phidiff_old = 0;
nn_delta = 0;
for nn = 1:1:num_t
    phiu = dt*ifft2( t(nn)*signMask.*exp(-t(nn).*pfs).*int);
    phi = phi + phiu;
    phidiff = sum( abs( phiu(:) - phiuo(:))); 
    %fprintf('%g ',phidiff)
    delta_phidiff = abs( phidiff_old - phidiff ); 
    if delta_phidiff < eps('single')
         nn_delta = nn_delta + 1;
         if nn > 1 && nn_delta > 10
            nn_break = nn;
         end
    end
    if nn > 1 && phidiff < eps('single')
        fprintf('\nAbort loop after %u iterations. Difference norm between subsequent updates: %g\n',nn,phidiff);
        nn_break = nn;        
    end    
    if (mod( nn, show_updates) == 0 ) || nn == 0 || nn_break > 0 || nn == num_t 
        subplot(2,1,1)
        imsc( phi )        
        title(sprintf('phi. iteration: %u', nn))
        colorbar
        axis equal tight
        
        subplot(2,1,2)
        imsc( phiu )
        title(sprintf('phi update. diff: %g, delta diff: %g', phidiff, delta_phidiff))
        colorbar
        axis equal tight  
        
        drawnow
        pause( 0.1 )
    end
    if nn == nn_break && nn > 1
        %keyboard
        break
    end
    phiuo = phiu;
    phidiff_old = phidiff;
    %phi(:,:,nn) = real(dt*ifft2( t(nn)*signMask.*exp(-t(nn).*pfs).*int));
end


figure

subplot(2,2,1)
imsc(phiRegOpt)
axis equal tight
title(sprintf('phi reg opt'))
colorbar
%xticks([]);yticks([])        

subplot(2,2,2)
imsc(phi)
axis equal tight
title(sprintf('phi'))
colorbar
%xticks([]);yticks([])        

h = @(x) 100*abs(x-phi0);

subplot(2,2,3)
imsc(h(phiRegOpt))
axis equal tight
title(sprintf('100*diff(phi reg opt,phi0)'))
colorbar
%xticks([]);yticks([])        

subplot(2,2,4)
imsc(h(phi))
axis equal tight
title(sprintf('100*diff(phi,phi0)'))
colorbar
%xticks([]);yticks([])        

drawnow

h = @(x) sum(abs(x(:)-phi0(:)));
fprintf('\nDifference norm: Standard: %g, Schwinger: %g\n',h(phiRegOpt),h(phi))

domain(phiRegOpt),domain(phi)
domain(signMask),domain(pfs),domain(imMask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signMask, phaseFilterSqrt, imMask] = SchwingerFilter(Method,imSize,EnergyDistancePixelsize,BinaryFilterThreshold)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    Method = 'tie';
end
if nargin < 3
    imSize = [1024 1024];
end


if nargin < 3
    EnergyDistancePixelsize = [20 0.945 .75e-6];
end
if nargin < 4
    BinaryFilterThreshold = 0.1;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
% wave length
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;

%% Fourier coordinates
xi  = FrequencyVector(imSize(2),'double',1);
eta = FrequencyVector(imSize(1),'double',1);
% 2D
[sinArg, sinxiquad]   = meshgrid(xi,eta);
% Function on 2D
sinArg = ArgPrefac*(sinArg.^2 + sinxiquad.^2)/2;

%% Filter 
switch lower(Method)
    case 'tie'
        signMask = 1/2*ones(imSize);
        phaseFilterSqrt = sqrt(sinArg);
        imMask = 1;
    case 'ctf'
        sinxiquad   = sin(sinArg);
        signMask    = 1/2*sign(sinxiquad);
        phaseFilterSqrt = sqrt(abs(sinxiquad));
        imMask     = 1;
    case {'qp','pctf','quasi','quasiparticle','quasiparticles'}
        sinxiquad   = sin(sinArg);
        signMask    = 1/2*sign(sinxiquad);
        phaseFilterSqrt = abs(sinxiquad);
        imMask     = ones(imSize);
        imMask( sinArg > pi/2  &  phaseFilterSqrt < BinaryFilterThreshold) = 0;
        phaseFilterSqrt = sqrt(phaseFilterSqrt);
    case {'qp2','qpmod','qpnew','quasinew','quasi2','pctf2','sinesupression'}
        %!!!!!!!!!!! TO CHECK !!!!!!!!!!!!!!!!!!!!!!!
        sinxiquad       = sin(sinArg);
        signMask        = 1/2*sign(sinxiquad);
        phaseFilterSqrt = abs(sinxiquad);
        mask            = sinArg > pi/2  &  phaseFilterSqrt < BinaryFilterThreshold;
        imMask         = ones(imSize);
        imMask( mask ) = bsxfun(@(a,b) a(b),2*signMask/(2*BinaryFilterThreshold),mask);
        phaseFilterSqrt = sqrt(phaseFilterSqrt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
