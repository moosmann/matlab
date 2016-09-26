function [out,g_ft,int,uncropped]=RecoGPU2(int,alpha,EnergyDistancePixelsize,evalTEIlo,evalTIEpnlo,evalCTF,BinaryFilter,Padding_FactorAndValue)
%TIE or CTF phase retrieval and additionaly Perturbatively evaluated
%Next-to-Leading Order (PNLO) correction to the TIE retrieved phase or extension of the
%CTF-applicability limit by means of a projection of the intput data onto
%the CTF model (projected CTF). Assuming pure-phase-contrast intensity maps with
%homogeneous absorption.
%
% int: flat-(and-dark-)field corrected radiograph.
% alpha: the regularization parameter is given as 10^-alpha with alpha in
% [9 18] for absolutly pure phase object (simulations). For slighty
% absorbing data alpha is usually in [3 5]. 
% EnergyDistancePixelsize: 3-vector, values in meter.
% BinaryFilter (for projeted CTF is): scalar with values in [0 1], for 0 no
% projection is computed.
% evalTIEpnlo: PNLO correction to TIE. 1- yes, 0 - no
% Padding_FactorAndValue: input is a cell array given in braces as {scalar
% string or scalar}. Padding value is either a number with which the data
% is padded or one of the predefined methods of Matlab's padarray function:
% 'circular', 'replicate' or 'symmetric'.
%
%[out,g_ft,int,uncropped]=RecoGPU2(int,alpha,EnergyDistancePixelsize,evalTEIlo,evalTIEpnlo,evalCTF,BinaryFilter,Padding_FactorAndValue)

%% Default arguments.
if nargin<2,alpha = 12;end
if nargin<3,EnergyDistancePixelsize = [25 1 0.36e-6];end % in units [keV m m]
if nargin<4,evalTEIlo = 1;end
if nargin<5,evalTIEpnlo = 0;end
if nargin<6,evalCTF = 1;end
if nargin<7,BinaryFilter = 0;end
if nargin<8,Padding_FactorAndValue = {1 'symmetric'};end

%% PARAMETERS.
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
% xi = matrix of identic row ranging accroding to the first entry of
% meshgrid, the rows are repeated according to the length of the vector
% of the second entry of meshgrid. eta analague to xi.
xi         = gpuArray(fftshift(xi));
eta        = gpuArray(fftshift(eta));
%xieta      = xi.*eta;
% Filter for CTF phase retrieval.
sineXiEta  = sin(1/prefactor*(xi.^2 + eta.^2)/2);
if BinaryFilter > 0 
    % Create binary filter with threshold 'BinaryFilter' for (1/prefactor*(xi.^2 + eta.^2)/2>pi/2)
    f = ones(dimx,dimy);
    f(gather(sineXiEta).^2<BinaryFilter) = 0;
    f(gather(xi.^2 + eta.^2)<pi*prefactor) = 1;
    f = gpuArray(f);
end
% Filter for linear TIE phase retrieval.
inv_sine   = 1./(2*sign(sineXiEta).*(abs(sineXiEta))+10^-alpha);
% Pad intensity to dimensions [dimx dimy].
int        = padarray(int,[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PaddingValue,'pre');
int        = padarray(int,[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PaddingValue,'post');
int        = gpuArray(int);
% Define intensity contrast g=I(x,y,z)/I(x,y,0)-1 and compute the Fourier
% transform of g. For pure phase objects. <g> states the conservation of
% flux <I>=1.  Substract mean is equivalent to the above definition of g in
% the case of pure pahse objects. To play safe the mean is substracted
% instead of 1. Since the mean value is not retrievable this doesn't change
% the phase retrieval. The mean value (or (1,1) component of g_ft) could
% be such large that it's beyond numerical precision.
g_ft       = int-sum(sum(int))/dimx/dimy;
g_ft       = g_ft-sum(sum(g_ft))/dimx/dimy;
g_ft       = fft2(g_ft);
%% CTF PHASE RETRIEVAL
% Divide FT[ctf] by the sine in Fourier space, invert
% the Fourier transform, take real part (just to play safe, substract
% mean (since it cannot be retrieved).
if evalCTF
    out.ctf       = real(ifft2(inv_sine.*g_ft));
    out.ctf       = gather(out.ctf-sum(sum(out.ctf))/dimx/dimy);
    if nargout > 2
        uncropped.ctf = out.ctf;
    end
    out.ctf       = out.ctf(xcut,ycut,:);
end
% Projected CTF
if BinaryFilter > 0,
    % Projected CTF: Before inverse FT apply binary filter to FT of CTF
    % phase.
    out.ctfProjected  = real(ifft2(f.*inv_sine.*g_ft));
    out.ctfProjected  = gather(out.ctfProjected-sum(sum(out.ctfProjected))/dimx/dimy);
    if nargout > 2
        uncropped.ctfProjected = out.ctfProjected;
    end
    out.ctfProjected  = out.ctfProjected(xcut,ycut,:);
end

%% TIE PHASE RETRIEVAL.
% Fourier transform of TIE retrieved phase. Here, denotetd out.tie to safe
% memory.
out.tie       = 1./(xi.^2 + eta.^2 + 10^-alpha).*g_ft;
% PNLO: Perturbatively evaluated Next-to-Leading Order correction.
if evalTIEpnlo > 0
    % Compute the correction to the laplacian of the phase, denoted as
    % out.tiePNLO to save variable.
    phi_dx2     = ifft2(eta.^2.*out.tie);
    phi_dy2     = ifft2( xi.^2.*out.tie);
    out.tiePNLO = -real( ...
        + ifft2(eta.*out.tie).*( ifft2( (    (eta.^2 + xi.^2).*eta).*out.tie)) ...
        + ifft2( xi.*out.tie).*( ifft2( (xi.*(eta.^2 + xi.^2)     ).*out.tie)) ...
        + phi_dx2.^2 ... 
        + phi_dx2.*phi_dy2 ... 
        + phi_dy2.^2 ...
        + ifft2(xi.*eta.*out.tie).^2);
    % Substract the mean of the Laplacian of the NLO.
    out.tiePNLO = out.tiePNLO - sum(sum(out.tiePNLO))/dimx/dimy;
    % Compute the PNLO of the phase by means of Inversion of the
    % Laplacian using inverse Fourier transforms.
    out.tiePNLO = ifft2(1./(xi.^2 + eta.^2 + 10^-alpha).*fft2(out.tiePNLO));
    % Take real part (just to play safe), rescale by prefactor, substract mean and crop back to
    % input dimensions.
    out.tiePNLO = prefactor*real(out.tiePNLO);
    out.tiePNLO = gather(out.tiePNLO-sum(sum(out.tiePNLO))/dimx/dimy);
    % Add uncropped map to output struct for debugging. Not necessary.
    if nargout > 2
        uncropped.tiePNLO = out.tiePNLO;
    end
    % The PNLO phase map.
    out.tiePNLO = out.tiePNLO(xcut,ycut,:);
end
%% Linear TIE: Leading order.
if evalTEIlo
    % Compute TIE phase by means of Fourier
    % inversion. Take real part (to play safe).
    out.tie      = prefactor*real(ifft2(out.tie));
    % Rescale by prefactor, substract mean and crop to
    % input dimensions.
    out.tie = gather(out.tie - sum(sum(out.tie))/dimx/dimy);
    % Add uncropped map to output struct for debugging. Not necessary.
    if nargout > 2
        uncropped.tie = out.tie;
    end
    % The LO-TIE phase map.
    out.tie = out.tie(xcut,ycut,:);
else
    out = rmfield(out,'tie');
end

% Add parameters to output struct 'out'.
out.edp                   = EnergyDistancePixelsize;
out.alpha                 = alpha;
out.BinaryFilterThreshold = BinaryFilter;
out.PaddedDimensionsXY    = [dimx dimy];
out.PaddingFactor         = PaddingFactor;
out.PaddingValue          = PaddingValue;
out.CroppedAreaXY         = {xcut; ycut};
