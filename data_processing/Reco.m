function [out, BinaryFilter, lffilt, uncropped] = Reco(int,alphaCTF_alphaTIE,EnergyDistancePixelsize,TIE_pCTF_CTF_PNLO,filterLowFreq,Padding_FactorMethod)
%TIE or CTF phase retrieval and additionaly Perturbatively evaluated
%Next-to-Leading Order (PNLO) correction to the TIE retrieved phase or extension of the
%CTF-applicability limit by means of a projection of the intput data onto
%the CTF model (projected CTF). Assuming pure-phase-contrast intensity maps with
%homogeneous absorption. Mostly, inevitable large-scale absorption
%contributions distort the retrieval of low-frequency components of the
%phase. This also affects the choice the value of the regularization
%parameter.
%
% int: 2D-matrix, flat-(and-dark-)field corrected radiograph.
% alphaCTF_alphaTIE: 1x2-vector, the regularization parameter is given as 10^-alpha with alpha in
% [9 18] for absolutly pure phase object (simulations). For slighty
% absorbing data a good choice of alpha is usually in the range [2 5], 2.5
% for CTF and 5 for TIE.
% EnergyDistancePixelsize: 1x3-vector, values in meter.
% TIE_pCTF_CTF_PNLO = [evalTIElo BinaryFilterThreshold evalCTF evalTIEpnlo]
%   evalTIElo: binary, compute linearized TIE phase. 1-yes, 0-no
%   BinaryFilterThreshold (for projeted CTF): scalar with values in [0 1] or
%      1x3-vector, for 0 (default) no projection is computed. If input is a
%      3-vector then a Gaussian-blurring is applied to the filter, 1st 
%      component is the threshold for the binary filter, 2nd component is the
%      support in pixels of the Gaussian filter (recommended: 3) and the 3rd 
%      component the sigma of which in pixels (recommended: 1)
%   evalCTF: binary, compute CTF phase. 1-yes, 0-no
%   evalTIEpnlo: binary, PNLO correction to TIE. 1-yes, 0-no
% Padding_FactorMethod: input is a cell array given in braces as {scalar
% string or scalar}. Padding value is either a number with which the data
% is padded or one of the predefined methods of Matlab's padarray function:
% 'circular', 'replicate' or 'symmetric'.
%
%[out BinaryFilter lffilt uncropped] = Reco(int,alphaCTF_alphaTIE,EnergyDistancePixelsize,TIE_pCTF_CTF_PNLO,filterLowFreq,Padding_FactorMethod)

%% Default arguments.
if nargin<2
    alphaCTF_alphaTIE = 2.5;
end
if nargin<3
    EnergyDistancePixelsize = [30 0.62 2.2e-6];
end % in units [keV m m]
if nargin<4
    TIE_pCTF_CTF_PNLO = [1 0 0 0];
end
if nargin < 5
    filterLowFreq = 30;
end
if nargin < 6
    Padding_FactorMethod = {1 'symmetric'};
end
%% PARAMETERS.
% Linearized TIE
evalTIElo = TIE_pCTF_CTF_PNLO(1);
% Projected CTF: Threshold for binary filter
if numel(TIE_pCTF_CTF_PNLO) > 1
    BinaryFilterThreshold = TIE_pCTF_CTF_PNLO(2);
else
    BinaryFilterThreshold = 0;
end
% CTF
if numel(TIE_pCTF_CTF_PNLO) > 2
    evalCTF = TIE_pCTF_CTF_PNLO(3);
else
    evalCTF = 0;
end
% Perturbatively evaluated next-to-leading order (NLO) correction to
% linearized TIE
if numel(TIE_pCTF_CTF_PNLO) > 3
    evalTIEpnlo = TIE_pCTF_CTF_PNLO(4);
else
    evalTIEpnlo = 0;
end
% Padding
PadFac  = Padding_FactorMethod{1};
PadMeth = Padding_FactorMethod{2};
% Experimetnal parameter
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
% Wave length.
lambda        = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
% Prefactor needed for TIE and CTF retrieval.
prefactor     = Pixelsize^2/(2*pi*lambda*Distance);
% Dimensions of input data array.
[dim1,dim2]   = size(int);
% Regularization parameters.
if size(alphaCTF_alphaTIE,2) == 1
    alphaCTF = alphaCTF_alphaTIE;
    alphaTIE = alphaCTF + log10(2*pi*lambda*Distance/Pixelsize^2);
else
    alphaCTF = alphaCTF_alphaTIE(1);
    alphaTIE = alphaCTF_alphaTIE(2);
end
% Get new dimension whichs are the next-power-of-2 times PadFac. The
% next-power-or-2 is done not to spoil computational power/precision because
% the fft of Matlab always pads to the next-power-of-2.
dimx          = PadFac*dim1;
dimy          = PadFac*dim2;
xcut          = 1+ceil((dimx-dim1)/2):floor((dimx+dim1)/2);
ycut          = 1+ceil((dimy-dim2)/2):floor((dimy+dim2)/2);
%% Filters.
% Fourier coordinates.
[xi,eta]   = meshgrid(-1/2:1/dimy:1/2-1/dimy,-1/2:1/dimx:1/2-1/dimx);
% xi = matrix of identic row ranging accroding to the first entry of
% meshgrid, the rows are repeated according to the length of the vector
% of the second entry of meshgrid. eta analague to xi.
xi  = fftshift(xi);
eta = fftshift(eta);
% Filter for CTF phase retrieval.
SineXiEta   = sin(1/prefactor*(xi.^2 + eta.^2)/2);
InverseSine = 1./(2*sign(SineXiEta).*(abs(SineXiEta))+10^-alphaCTF);
%% Padding and FT of data.
% Pad intensity to dimensions [dimx dimy].
int        = padarray(int,[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PadMeth,'pre');
int        = padarray(int,[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PadMeth,'post');
if nargout > 3
    uncropped.intensity = int;
end
% Define intensity contrast g=I(x,y,z)/I(x,y,0)-1 and compute the Fourier
% transform of g. For pure phase objects. <g> states the conservation of
% flux <I>=1.  Substract mean is equivalent to the above definition of g in
% the case of pure pahse objects. To play safe the mean is substracted
% instead of 1. Since the mean value is not retrievable this doesn't change
% the phase retrieval. The mean value (or (1,1) component of g_ft) could
% be such large that it's beyond numerical precision.
%% !! FT of 'int' has same variable name to save memory.
int       = int-mean(int(:));
int       = fft2(int);
if nargout > 3
    uncropped.intensityFT = int;
end
%% Filter low requencies
if filterLowFreq > 0
    PadFacX       = dimx/dim1;
    PadFacY       = dimy/dim2;
    lffilt = (1-exp(-((dimy*xi/PadFacY/filterLowFreq).^2+(dimx*eta/PadFacX/filterLowFreq).^2)/2));
    int = lffilt.*int;
end
%% CTF PHASE RETRIEVAL
% Divide FT[ctf] by the sine in Fourier space, invert
% the Fourier transform, take real part (just to play safe, substract
% mean (since it cannot be retrieved).
if evalCTF
    out.ctf       = real(ifft2(InverseSine.*int));
    out.ctf       = out.ctf-mean(out.ctf(:));
    if nargout > 3
        uncropped.ctf = out.ctf;
    end
    out.ctf       = out.ctf(xcut,ycut);
end
%% Projected CTF
if BinaryFilterThreshold(1) > 0,
    % Create binary filter with threshold 'BinaryFilterThreshold(1)' for (1/prefactor*(xi.^2 + eta.^2)/2>pi/2)
    BinaryFilter = ones(dimx,dimy);
    BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
    if size(BinaryFilterThreshold,2) > 1
        hsize = BinaryFilterThreshold(2);
        sigma = BinaryFilterThreshold(3);
        BinaryFilter = imfilter(BinaryFilter,fspecial('gaussian',[hsize hsize],sigma));
        BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)/10) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
    end
    % Projected CTF: Before inverse FT apply binary filter to FT of CTF
    % phase.
    out.ctfProjected  = real(ifft2(BinaryFilter.*InverseSine.*int));
    out.ctfProjected  = out.ctfProjected-mean(out.ctfProjected(:));
    if nargout > 3
        uncropped.ctfProjected = out.ctfProjected;
    end
    out.ctfProjected  = out.ctfProjected(xcut,ycut);
    % If BinaryFilter is wanted for output uncomment next line.
    %BinaryFilter  = fftshift(BinaryFilter);
end
%% TIE PNLO (Perturbatively evaluated Next-to-Leading Order)
% Fourier transform of TIE retrieved phase. Here, denotetd out.tieLO to safe
% memory.
if evalTIElo || evalTIEpnlo
    out.tieLO       = 1./(xi.^2 + eta.^2 + 10^-alphaTIE).*int;
end
% PNLO: Perturbatively evaluated Next-to-Leading Order correction.
if evalTIEpnlo > 0
    % Compute the correction to the laplacian of the phase, denoted as
    % out.tieLOPNLO to save variable.
    phi_dx2     = ifft2(eta.^2.*out.tieLO);
    phi_dy2     = ifft2( xi.^2.*out.tieLO);
    out.tiePNLO = -real( ...
        + ifft2(eta.*out.tieLO).*( ifft2( (    (eta.^2 + xi.^2).*eta).*out.tieLO)) ...
        + ifft2( xi.*out.tieLO).*( ifft2( (xi.*(eta.^2 + xi.^2)     ).*out.tieLO)) ...
        + phi_dx2.^2 ...
        + phi_dx2.*phi_dy2 ...
        + phi_dy2.^2 ...
        + ifft2(xi.*eta.*out.tieLO).^2);
    % Substract the mean of the Laplacian of the NLO.
    out.tiePNLO = out.tiePNLO - mean(out.tiePNLO(:));
    % Compute the PNLO of the phase by means of Inversion of the
    % Laplacian using inverse Fourier transforms.
    out.tiePNLO = ifft2(1./(xi.^2 + eta.^2 + 10^-alphaTIE).*fft2(out.tiePNLO));
    % Take real part (just to play safe), rescale by prefactor, substract mean and crop back to
    % input dimensions.
    out.tiePNLO = prefactor*real(out.tiePNLO);
    out.tiePNLO = out.tiePNLO-mean(out.tiePNLO(:));
    % Add uncropped map to output struct for debugging. Not necessary.
    if nargout > 3
        uncropped.tiePNLO = out.tiePNLO;
    end
    % The PNLO phase map.
    out.tiePNLO = out.tiePNLO(xcut,ycut);
    if ~evalTIElo
        out = rmfield(out,'tieLO');
    end
end
%% Linear TIE (Leading order).
if evalTIElo
    % Compute TIE phase by means of Fourier
    % inversion. Take real part (to play safe).
    out.tieLO      = real(ifft2(out.tieLO));
    % Rescale by prefactor, substract mean and crop to
    % input dimensions.
    out.tieLO = out.tieLO - mean(out.tieLO(:));
    out.tieLO = prefactor*out.tieLO;
    out.tieLO = out.tieLO - mean(out.tieLO(:));
    % Add uncropped map to output struct for debugging. Not necessary.
    if nargout > 3
        uncropped.tie = out.tieLO;
    end
    % The LO-TIE phase map.
    out.tieLO = out.tieLO(xcut,ycut);
end
if nargout > 1
    BinaryFilter = fftshift(BinaryFilter);
end
%% Add parameters to output struct 'out'.
out.edp                   = EnergyDistancePixelsize;
out.alpha                 = alphaCTF_alphaTIE;
out.BinaryFilterThreshold = BinaryFilterThreshold;
out.PaddedDimensionsXY    = [dimx dimy];
out.PaddingFactor         = PadFac;
out.PaddingMethod          = PadMeth;
out.CroppedAreaXY         = {xcut; ycut};
