function Schwinger2(im,Method,EnergyDistancePixelsize,BinaryFilterThreshold)
% Another try to solve the regularization problem using the Schwinger
% inegral to circumvent the singularity of 1/x^2 at x = 0.
%
% Written by Julian Moosmann, first version: 2013-12-19, last version: 2013-12-19
%
% Schwinger2()

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    Method = 'qp';    
end
if nargin < 3
    EnergyDistancePixelsize = [30 0.7 1e-6];
end
if nargin < 4
    BinaryFilterThreshold = 0.05;
end
if nargin < 1
    phi0 = 0.1*normat(double(imread('~/barbara.png')));
    im = Propagation2(phi0,0,EnergyDistancePixelsize,0,1);
    phi0 = SubtractMean(phi0);
end
if nargin < 5
    loopOverRegPar = 1;
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FT of intensity
im = fft2(im-1);


% reference retrieval
if loopOverRegPar(1)
    h = @(x) x(:);
    regPar = 2:0.2:8;
    for mm = numel(regPar):-1:1
        phiReg(:,:,mm) = ifft2( PhaseFilter(Method,size(im),EnergyDistancePixelsize,regPar(mm),0.1,'double') .* im );
        er(mm) = sum( abs( h( phiReg(:,:,mm) )  - phi0(:) ) );
        er2(mm) = sum( abs( h( normat(phiReg(:,:,mm) ) ) - h( normat(phi0) ) ) );
    end
    [erMinVal, erMinPos] = min(er(:));
    fprintf('\n Minimum of error map of standard regularized retrieval at index %u, regPar %g, (and error value %g)',erMinPos,regPar(erMinPos),erMinVal)
    [er2MinVal, er2MinPos] = min(er2(:));
    fprintf('\n Minimum of error map of renormed standard regularized retrieval at index %u, regPar %g, (and error value %g)\n',er2MinPos,regPar(er2MinPos),er2MinVal)
    phiRegOpt = phiReg(:,:,er2MinPos);
else
    phiRegOpt = ifft2( PhaseFilter(Method,size(im),EnergyDistancePixelsize,4.4,BinaryFilterThreshold,'double') .* im );    
end

% create filters and masks for Schwinger integration
[signMask, pfs, intMask] = SchwingerFilter(Method,size(im),EnergyDistancePixelsize,BinaryFilterThreshold);

% mask for QP
im = intMask .* im;

% determine tMax where the exponential weight drops below numerical
% precision, the smallest values of  pfs=sqrt(abs(sin(pi lambda xi^2 /
% dx^2))) is close the value of pfs(2,1)=pfs(1,2), pfs(1,1) = 0. There
% are/could be smaller values close to higher order of zeros, but the order
% of magnitude should be the same
tMax = -log(eps)/pfs(2,1);
dt = 0.2;0.25;%tMax/1000;
t = dt:dt:tMax;
% integrate using a loop
phi = 0;
phiuo = 0;
%phi = zeros([size(im) numel(t)]);
for nn = 1:1:numel(t)
    phiu = dt*ifft2( t(nn)*signMask.*exp(-t(nn).*pfs).*im);
    phi = phi + phiu;
    phidiff = sum( abs( phiu(:) - phiuo(:))); 
    %fprintf('%g ',phidiff)
    if phidiff < 10^-6
        fprintf('\nAbort loop after %u iterations. Difference norm between subsequent updates: %g\n',nn,phidiff);
        break
    end
    phiuo = phiu;
    %phi(:,:,nn) = real(dt*ifft2( t(nn)*signMask.*exp(-t(nn).*pfs).*im));
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
domain(signMask),domain(pfs),domain(intMask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signMask, phaseFilterSqrt, intMask] = SchwingerFilter(Method,imSize,EnergyDistancePixelsize,BinaryFilterThreshold)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    Method = 'tie';
end
if nargin < 2
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
        intMask = 1;
    case 'ctf'
        sinxiquad   = sin(sinArg);
        signMask    = 1/2*sign(sinxiquad);
        phaseFilterSqrt = sqrt(abs(sinxiquad));
        intMask     = 1;
    case {'qp','pctf','quasi','quasiparticle','quasiparticles'}
        sinxiquad   = sin(sinArg);
        signMask    = 1/2*sign(sinxiquad);
        phaseFilterSqrt = abs(sinxiquad);
        intMask     = ones(imSize);
        intMask( sinArg > pi/2  &  phaseFilterSqrt < BinaryFilterThreshold) = 0;
        phaseFilterSqrt = sqrt(phaseFilterSqrt);
    case {'qp2','qpmod','qpnew','quasinew','quasi2','pctf2','sinesupression'}
        %!!!!!!!!!!! TO CHECK !!!!!!!!!!!!!!!!!!!!!!!
        sinxiquad       = sin(sinArg);
        signMask        = 1/2*sign(sinxiquad);
        phaseFilterSqrt = abs(sinxiquad);
        mask            = sinArg > pi/2  &  phaseFilterSqrt < BinaryFilterThreshold;
        intMask         = ones(imSize);
        intMask( mask ) = bsxfun(@(a,b) a(b),2*signMask/(2*BinaryFilterThreshold),mask);
        phaseFilterSqrt = sqrt(phaseFilterSqrt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
