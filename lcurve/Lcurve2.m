function [out NormOfSolution NormOfResidualOfSolution] = Lcurve2(int,AlphaIntervall,NumOfPts,EnergyDistancePixelsize,doLogOfNorm,Padding_FactorAndValue)

% Compute Lcurve to determine the regularization parameter for the phase
% retrieval problem.

%% Defaults arguments.
if nargin<2,AlphaIntervall=[0 5];end
if nargin<3,NumOfPts=11;end
if nargin<4,EnergyDistancePixelsize=[1 1 1];end
if nargin<5,doLogOfNorm=1;end;
if nargin<6,Padding_FactorAndValue={1 'symmetric'};end

%% PARAMETERS.
AlphaMin = AlphaIntervall(1);
AlphaMax = AlphaIntervall(2);
PaddingFactor = Padding_FactorAndValue{1};
PaddingValue  = Padding_FactorAndValue{2};
% Wave length.
lambda        = 6.62606896e-34*299792458/(EnergyDistancePixelsize(1)*1.60217733e-16);
% Prefactor needed for TIE and CTF retrieval.
prefactor     = EnergyDistancePixelsize(3)^2/(2*pi*lambda*EnergyDistancePixelsize(2));
% Dimensions of input data array.
[dim1,dim2]   = size(int);
% Get new dimension whichs are the next-power-of-2 times PaddingFactor. The
% next-power-or-2 is done not to spoil computational power/precision because
% the fft of Matlab always pads to the next-power-of-2.
dimx          = PaddingFactor*2^nextpow2(dim1);
dimy          = PaddingFactor*2^nextpow2(dim2);
xcut          = 1+ceil((dimx-dim1)/2):floor((dimx+dim1)/2);
ycut          = 1+ceil((dimy-dim2)/2):floor((dimy+dim2)/2);
%% PROGRAMME.
% Fourier coordinates.
[xi,eta]   = meshgrid(-1/2:1/dimy:1/2-1/dimy,-1/2:1/dimx:1/2-1/dimx);
xi         = fftshift(xi);
eta        = fftshift(eta);
lap        = xi.^2 + eta.^2;
% Pad intensity to dimensions [dimx dimy].
int        = padarray(int,[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PaddingValue,'pre');
int        = padarray(int,[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PaddingValue,'post');
% Define intensity contrast g=I(x,y,z)/I(x,y,0)-1 and compute the Fourier
% transform of g.
int        = int-mean(int(:));
int        = int-mean(int(:));
gFT        = fft2(int);
%% L-curve: Loop over regularization parameter values.
tic;
for nn=NumOfPts:-1:1
    % Regularization parameter m^2=10^-Alpha
    Alpha(nn) = AlphaMin + (AlphaMax-AlphaMin)*(nn-1)/(NumOfPts-1);
    % Linear TIE: Leading order.
    % Compute TIE phase by means of Fourier
    % inversion. Take real part (to play safe) and rescale.
    phaseFT     = prefactor*1./(lap + 10^-Alpha(nn)).*gFT;
    out(nn).tie = real(ifft2(phaseFT));
    %fprintf('Norm: %10g, Mean: %10g\n',norm(out(nn).tie(:)),mean(out(nn).tie(:)));
    out(nn).tie = out(nn).tie - mean(out(nn).tie(:));   
    %fprintf('Norm: %10g, Mean: %10g\n',norm(out(nn).tie(:)),mean(out(nn).tie(:)));
    % Norm of the solution and norm of the residual of the Solution.
    phaseFT     = 1/prefactor*real(ifft2(lap.*phaseFT));
    %out(nn).IntFromData = phaseFT(xcut,ycut);
    phaseFT     = phaseFT - int;
    switch doLogOfNorm
        case {1,'yes'}
            NormOfSolution(nn)           = log10(norm(out(nn).tie(:)));
            NormOfResidualOfSolution(nn) = log10(norm(phaseFT(:)));
        case {0,'no'}
            NormOfSolution(nn)           = norm(out(nn).tie(:));
            NormOfResidualOfSolution(nn) = norm(phaseFT(:));
    end    
    % Crop phase map back to original size.
    out(nn).tie = out(nn).tie(xcut,ycut);
    % Add parameters to output struct 'out'.
    out(nn).edp                      = EnergyDistancePixelsize;
    out(nn).Alpha                    = Alpha(nn);
    out(nn).PaddedDimensionsXY       = [dimx dimy];
    out(nn).PaddingFactor            = PaddingFactor;
    out(nn).PaddingValue             = PaddingValue;
    out(nn).CroppedAreaXY            = {xcut; ycut};
    out(nn).NormOfSolution           = NormOfSolution(nn);
    out(nn).NormOfResidualOfSolution = NormOfResidualOfSolution(nn);
    fprintf('Index: %2u, Alpha = %4f, Norm of the solution: %8.3g, Norm of the residual of the solution: %8.3g\n', ... 
        nn,Alpha(nn),NormOfSolution(nn),NormOfResidualOfSolution(nn));
end
tReco = toc;
%% Compute curvature of L-curve.
NormSol_D1 = diff(abs(NormOfSolution),1);
NormSol_D1 = interp1(1.5:NumOfPts,NormSol_D1,2:NumOfPts-1);
NormSol_D2 = diff(abs(NormOfSolution),2);
ResSol_D1  = diff(abs(NormOfResidualOfSolution),1);
ResSol_D1  = interp1(1.5:NumOfPts,ResSol_D1,2:NumOfPts-1);
ResSol_D2  = diff(abs(NormOfResidualOfSolution),2);
curvature  = abs((ResSol_D1.*NormSol_D2 - ResSol_D2.*NormSol_D1)./(NormSol_D1.^2 + ResSol_D2.^2).^1.5);
[curvatureMax nn_curvatureMax] = max(curvature(2:end-1));
nn_curvatureMax = nn_curvatureMax + 2;

%% Print info and plots.
fprintf('Linearized TIE phase retrieval for %u regularization parameter values done in %gs.\n',NumOfPts,tReco)
switch doLogOfNorm
    case {1,'yes'}
        fprintf('Took log of norm and residual of the solution.\n')
    case {0,'no'}
        fprintf('Did NOT took log of norm and residual of the solution\n')
end
fprintf('At index %u and Alpha = %f the curvature has maximum of value %g.\n', ...
    nn_curvatureMax,AlphaMin+(AlphaMax-AlphaMin)*(nn_curvatureMax-1)/(NumOfPts-1),curvatureMax);
% Plots.
figure('Name','L-curve of linearized TIE phase retrieval (leading order)'),plot(NormOfResidualOfSolution,NormOfSolution,'.-'), ...
for nn=1:NumOfPts
    text(NormOfResidualOfSolution(nn),NormOfSolution(nn),sprintf('  ( %u , %.3g , %.3g )',nn,Alpha(nn),10^-Alpha(nn)));
end
figure('Name','Curvature of the L-curve vs regularization parameter'),plot(Alpha(2:NumOfPts-1),curvature);

end

