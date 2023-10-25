function [out,g_ft,uncropped]=RecoExplicitTIE(int,alpha,EnergyDistancePixelsize,BinaryFilter,PerturbativeCorrection,Padding_FactorAndValue)

% [out,g_ft,ctf_uncropped] =
% Reco(int,alpha,EnergyDistancePixelsize,BinaryFilter,PerturbativeCorrection,Padding_FactorAndValue);
                                                 
%% Default arguments.
if nargin<2,alpha=12;end
if nargin<3,EnergyDistancePixelsize=[25 1 0.36e-6];end
if nargin<4,BinaryFilter=0;end
if nargin<5,PerturbativeCorrection=0;end
if nargin<6,Padding_FactorAndValue={1 'symmetric'};end

%% PARAMETERS.
PaddingFactor = Padding_FactorAndValue{1};
PaddingValue  = Padding_FactorAndValue{2};
prefactor     = EnergyDistancePixelsize(3)^2/(2*pi*EnergyConverter(EnergyDistancePixelsize(1))*EnergyDistancePixelsize(2));
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
% Filters.
[xi,eta]   = meshgrid(-1/2:1/dimy:1/2-1/dimy,-1/2:1/dimx:1/2-1/dimx);
% xi = matrix of identic row ranging accroding to the first entry of
% meshgrid, the rows are repeated according to the length of the vector
% of the second entry of meshgrid. eta analague to xi.
xi         = fftshift(xi);
eta        = fftshift(eta);
%xieta      = xi.*eta;
sineXiEta  = sin(1/prefactor*(xi.^2 + eta.^2)/2);
inv_sine   = 1./(2*sign(sineXiEta).*(abs(sineXiEta))+10^-alpha);
% Pad intensity to dimensions [dimx dimy]: to avoid changing int. Define intensity contrast
% g=I(x,y,z)/I(x,y,0)-1. Compute the Fourier transform of g. For pure phase
% objects. <g> states the conservation of flux <I>=1.
int        = padarray(int,[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PaddingValue,'pre');
int        = padarray(int,[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PaddingValue,'post');
g_ft       = int-mean(int(:));
g_ft       = g_ft-mean(g_ft(:));
g_ft       = fft2(g_ft);

%% CTF PHASE RETRIEVAL
% Divide FT[ctf] by the sine in Fourier space, invert
% the Fourier transform, take real part (just to play safe, substract
% mean (since it cannot be retrieved).
out.ctf       = real(ifft2(inv_sine.*g_ft));
out.ctf       = out.ctf-mean(out.ctf(:));
if nargout == 3
    uncropped.ctf = out.ctf;
end
out.ctf       = out.ctf(xcut,ycut,:);
% Projected CTF
if BinaryFilter>0
    % Create binary filter with threshold BinaryFilter.
    f = ones(dimx,dimy);
    f( (sineXiEta.^2<BinaryFilter) & (1/prefactor*(xi.^2 + eta.^2)/2>pi/2) ) = 0;
    % Projected CTF: Before inverse FT apply binary filter to FT of CTF
    % phase.
    out.ctfProjected  = real(ifft2(f.*inv_sine.*g_ft));
    out.ctfProjected  = out.ctfProjected-mean(out.ctfProjected(:));
    if nargout == 3
        uncropped.ctfProjected = out.ctfProjected;
    end
    out.ctfProjected  = out.ctfProjected(xcut,ycut,:);
end

%% TIE PHASE RETRIEVAL.
% Fourier transform of TIE retrieved phase.
phi_ft       = 1./(xi.^2 + eta.^2 + 10^-alpha).*g_ft;
% PNLO: Perturbatively evaluated Next-to-Leading Order correction.
if PerturbativeCorrection > 0
    phi_dx1     = ifft2(eta   .*phi_ft);
    phi_dx2     = ifft2(eta.^2.*phi_ft);
    phi_dx3     = ifft2(eta.^3.*phi_ft);
    phi_dy1     = ifft2( xi   .*phi_ft);
    phi_dy2     = ifft2( xi.^2.*phi_ft);
    phi_dy3     = ifft2( xi.^3.*phi_ft);
    phi_dx1dy1  = ifft2(xi.*eta.*phi_ft);
    phi_dx1dy2  = ifft2( xi.^2.*eta.*phi_ft);
    phi_dx2dy1  = ifft2( xi.*eta.^2.*phi_ft);
    % Compute the correction to the laplacian of the phase, denote as
    % out.tiePNLO to save variable.
    out.tiePNLO = -real( ...
        + phi_dx1.*(phi_dx3    + phi_dx1dy2) ...
        + phi_dy1.*(phi_dx2dy1 + phi_dy3) ...
        + phi_dx2.^2 ... 
        + phi_dx2.*phi_dy2 ... 
        + phi_dx1dy1.^2 ...
        + phi_dy2.^2);
    % Substract the mean of the Laplacian of the correction of phase.
    out.tiePNLO = out.tiePNLO - mean(out.tiePNLO(:));
    % Subsract the correction to the phase by means of Inversion of the
    % Laplacian using Fourier transforms.
    out.tiePNLO = ifft2(1./(xi.^2 + eta.^2 + 10^-alpha).*fft2(out.tiePNLO));
    % Take real part (just to play safe), rescale, substract mean and crop back to
    % input dimensions.
    out.tiePNLO = prefactor*real(out.tiePNLO);
    out.tiePNLO = out.tiePNLO-mean(out.tiePNLO(:));
    if nargout == 3
        uncropped.tiePNLO = out.tiePNLO;
    end
    out.tiePNLO = out.tiePNLO(xcut,ycut,:);
end

% Linear TIE: Leading order.
% Compute TIE phase by means of Fourier inversion. Take real part (just
% to play safe).
out.tie      = real(ifft2(phi_ft));
% Rescale, substract mean and crop back to
% input dimensions.
out.tie = out.tie - mean(out.tie(:));
out.tie = prefactor*out.tie;
out.tie = out.tie - mean(out.tie(:));
if nargout == 3
    uncropped.tie = out.tie;
end
out.tie = out.tie(xcut,ycut,:);

% Add (input and computational) parameters to output struct 'out'.
out.edp                   = EnergyDistancePixelsize;
out.alpha                 = alpha;
out.BinaryFilterThreshold = BinaryFilter;
out.PaddedDimensionsXY    = [dimx dimy];
out.PaddingFactor         = PaddingFactor;
out.PaddingValue          = PaddingValue;
out.CroppedAreaXY         = {xcut; ycut};
