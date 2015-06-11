function [out,uncropped]=RecoGPU(Intensity,alpha,EnergyDistancePixelsize,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,Padding_FactorAndValue)
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
% BinaryFilterThreshold (for projeted CTF): scalar with values in [0 1],
% for 0 (default) no projection is computed.
% evalTIEpnlo: PNLO correction to TIE. 1- yes, 0 - no
% Padding_FactorAndValue: input is a cell array given in braces as {scalar
% string or scalar}. Padding value is either a number with which the data
% is padded or one of the predefined methods of Matlab's padarray function:
% 'circular', 'replicate' or 'symmetric'.
%
%[out,uncropped]=RecoGPU(Intensity,alpha,EnergyDistancePixelsize,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,Padding_FactorAndValue)

%% Default arguments.
if nargin<2,alpha = 12;end
if nargin<3,EnergyDistancePixelsize = [25 1 0.36e-6];end % in units [keV m m]
if nargin<4,evalTIElo = 1;end
if nargin<5,evalTIEpnlo = 0;end
if nargin<6,evalCTF = 1;end
if nargin<7,BinaryFilterThreshold = 0;end
if nargin<8,Padding_FactorAndValue = {1 'symmetric'};end

%% PARAMETERS.
PaddingFactor = Padding_FactorAndValue{1};
PaddingValue  = Padding_FactorAndValue{2};
% Wave length.
lambda        = 6.62606896e-34*299792458/(EnergyDistancePixelsize(1)*1.60217733e-16);
% Prefactor needed for TIE and CTF retrieval.
prefactor     = EnergyDistancePixelsize(3)^2/(2*pi*lambda*EnergyDistancePixelsize(2));
% Dimensions of input data array.
[dim1,dim2]   = size(Intensity);
% Get new dimension whichs are the next-power-of-2 times PaddingFactor. The
% next-power-or-2 is done not to spoil computational power/precision because
% the fft of Matlab always pads to the next-power-of-2.
dimx          = PaddingFactor*2^nextpow2(dim1);
dimy          = PaddingFactor*2^nextpow2(dim2);
xcut          = 1+ceil((dimx-dim1)/2):floor((dimx+dim1)/2);
ycut          = 1+ceil((dimy-dim2)/2):floor((dimy+dim2)/2);
%% Filters.
% Fourier coordinates.
[xi,eta]   = meshgrid(-1/2:1/dimy:1/2-1/dimy,-1/2:1/dimx:1/2-1/dimx);
% xi = matrix of identic row ranging accroding to the first entry of
% meshgrid, the rows are repeated according to the length of the vector
% of the second entry of meshgrid. eta analague to xi.
xi         = fftshift(xi);
eta        = fftshift(eta);
%xieta      = xi.*eta;
% Filter for CTF phase retrieval.
SineXiEta  = sin(1/prefactor*(xi.^2 + eta.^2)/2);
if BinaryFilterThreshold(1) > 0
    % Create binary filter with threshold 'BinaryFilterThreshold' for (1/prefactor*(xi.^2 + eta.^2)/2>pi/2)
    BinaryFilter = ones(dimx,dimy);
    %BinaryFilter(SineXiEta.^2<BinaryFilterThreshold) = 0;
    %BinaryFilter((xi.^2 + eta.^2)<pi*prefactor) = 1;
    BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
    if size(BinaryFilterThreshold,2) > 1
        hsize = BinaryFilterThreshold(2);
        sigma = BinaryFilterThreshold(3);
        BinaryFilter = imfilter(BinaryFilter,fspecial('gaussian',[hsize hsize],sigma));
        BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)/10) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
    end
    BinaryFilter = gpuArray(BinaryFilter);
end
% Filter for linear TIE phase retrieval.
InverseSine   = 1./(2*sign(SineXiEta).*(abs(SineXiEta))+10^-alpha);
InverseSine   = gpuArray(InverseSine);
xi         = gpuArray(xi);
eta        = gpuArray(eta);
% Pad intensity to dimensions [dimx dimy].
Intensity        = padarray(Intensity,[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PaddingValue,'pre');
Intensity        = padarray(Intensity,[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PaddingValue,'post');
if nargout > 1
    uncropped.intensity = Intensity;
end
% Define intensity contrast g=I(x,y,z)/I(x,y,0)-1 and compute the Fourier
% transform of g. For pure phase objects. <g> states the conservation of
% flux <I>=1.  Substract mean is equivalent to the above definition of g in
% the case of pure pahse objects. To play safe the mean is substracted
% instead of 1. Since the mean value is not retrievable this doesn't change
% the phase retrieval. The mean value (or (1,1) component of g_ft) could
% be such large that it's beyond numerical precision.
%% !! FT of 'Intensity' has same variable name to save memory.
Intensity       = Intensity-mean(Intensity(:));
Intensity       = gpuArray(Intensity);
Intensity       = fft2(Intensity);
if nargout > 1
    uncropped.intensityFT = Intensity;
end
%% CTF PHASE RETRIEVAL
% Divide FT[ctf] by the sine in Fourier space, invert
% the Fourier transform, take real part (just to play safe, substract
% mean (since it cannot be retrieved).
if evalCTF
    out.ctf       = gather(real(ifft2(InverseSine.*Intensity)));
    out.ctf       = out.ctf-mean(out.ctf(:));
    if nargout > 1
        uncropped.ctf = out.ctf;
    end
    out.ctf       = out.ctf(xcut,ycut);
end
% Projected CTF
if BinaryFilterThreshold > 0,
    % Projected CTF: Before inverse FT apply binary filter to FT of CTF
    % phase.
    out.ctfProjected  = gather(real(ifft2(BinaryFilter.*InverseSine.*Intensity)));
    out.ctfProjected  = out.ctfProjected-mean(out.ctfProjected(:));
    if nargout > 1
        uncropped.ctfProjected = out.ctfProjected;
    end
    out.ctfProjected  = out.ctfProjected(xcut,ycut);
    % If BinaryFilter is needed for output uncomment next line.
    %BinaryFilter  = fftshift(gather(BinaryFilter));
end
%% TIE PNLO (Perturbatively evaluated Next-to-Leading Order)
% Fourier transform of TIE retrieved phase. Here, denotetd out.tieLO to safe
% memory.
if evalTIElo || evalTIEpnlo
    out.tieLO       = 1./(xi.^2 + eta.^2 + 10^-alpha).*Intensity;
end
% PNLO: Perturbatively evaluated Next-to-Leading Order correction.
if evalTIEpnlo > 0
    % Compute the correction to the laplacian of the phase, denoted as
    % out.tiePNLO to save variable.
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
    %out.tiePNLO = out.tiePNLO - mean(out.tiePNLO(:));
    % Compute the PNLO of the phase by means of Inversion of the
    % Laplacian using inverse Fourier transforms.
    out.tiePNLO = ifft2(1./(xi.^2 + eta.^2 + 10^-alpha).*fft2(out.tiePNLO));
    % Take real part (just to play safe), rescale by prefactor, substract mean and crop back to
    % input dimensions.
    out.tiePNLO = gather(prefactor*real(out.tiePNLO));
    out.tiePNLO = out.tiePNLO-mean(out.tiePNLO(:));
    % Add uncropped map to output struct for debugging. Not necessary.
    if nargout > 1
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
    out.tieLO      = gather(real(ifft2(out.tieLO)));
    % Rescale by prefactor, substract mean and crop to
    % input dimensions.
    out.tieLO = out.tieLO - mean(out.tieLO(:));
    out.tieLO = prefactor*out.tieLO;
    out.tieLO = out.tieLO - mean(out.tieLO(:));
    % Add uncropped map to output struct for debugging. Not necessary.
    if nargout > 1
        uncropped.tie = out.tieLO;
    end
    % The LO-TIE phase map.
    out.tieLO = out.tieLO(xcut,ycut);
end
%% Add parameters to output struct 'out'.
out.edp                   = EnergyDistancePixelsize;
out.alpha                 = alpha;
out.BinaryFilterThreshold = BinaryFilterThreshold;
out.PaddedDimensionsXY    = [dimx dimy];
out.PaddingFactor         = PaddingFactor;
out.PaddingValue          = PaddingValue;
out.CroppedAreaXY         = {xcut; ycut};
