function [sinotie sinopctf] = RecoLoopThis2(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,EnergyPixelsizeDistance,InputPath,FileNamePrefix,OutputPath,doSino,IndexArray,DataSetName,ycrop,ProjDecr)
% Script for phase retrieval. Script loops over all folders contained in
% the CURRENT 'parent folder' and starting with 'FileNamePrefix'.

%% Default parameters.
if nargin < 1
    % Scalar: regularization parameter needed to fix the singularity encountered
    % when inverting the laplacian (TIE) or sine function (CTF).
    alphaCTF_alphaTIE  = 2.5;
end
if nargin < 2
    % Scalar: evaluate leading-order (LO) transport-of-intensity expression (TIE):
    % TIElo.
    evalTIElo    = 1;
end
if nargin < 3
    % Scalar: evaluate perturbatively next-to-leading-order (PNLO) TIE expression: TIEpnlo
    evalTIEpnlo  = 0;
end
if nargin < 4
    % Scalar: evaluate contrast-transfer-function (CTF) expression.
    evalCTF      = 0.1;
end
if nargin < 5
    % Scalar, or 1x3-Vector [Threshold Hsize Sigma]: If Threshold > 0:
    % evaluate projected CTF (PCTF) where the value is the threshold
    % for the projection filter. Optionally two more values can be passed
    % defining the size in pixels 'Hsize' and the standard deviation
    % 'Sigma' for a Gaussian blur filter which will be applied to the CTF
    % filter to smooth the edges (0/1-jumps) of which.
    BinaryFilterThreshold = 0;
end
if nargin < 6
    % 1x3-Vector, or cell of 1x3-Vectors: 1x3-Vector = [Energy Distance
    % Pixelsize] in metre. Can be cell array to specify for each data set
    % its own parameters. E.g.:
    EnergyPixelsizeDistance = [30.00 2.2e-6 0.620];
end
if nargin < 7
    % String: trailing seperator ('/') not needed. E.g.
    InputPath = '/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266/data/wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV/';
end
if nargin < 8
    % Prefix of filenamse under 'InputPath'
    FileNamePrefix = 'proj';
end
if nargin < 9
    OutputPath = '/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266/phase/wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV/';
end
if nargin < 10
    % Make and save sinograms of the center of projections for input to
    % inverse Radon transfomr (iradon). doSino > 1: slice number 'doSino'
    % is taken for the sinogram
    doSino = 1;
end
 % For HDF files.
if nargin < 11
     IndexArray  = {[1  1],[1  1],[1024  2048]};
end
if nargin < 12
    DataSetName = '/entry1/data/data';
end
if nargin < 13 
    % Vertical cropping.
    ycrop = 1024;%1024;
    xcrop = 1008;%128;
end
if nargin < 14
    ProjDecr = 0;
end
Padding_FactorAndValue = {1 'symmetric'};
StartProj = 50;
NumProj   = 1200;
EndProj   = StartProj + NumProj;
HotPixThrDark = 0.02;
HotPixThrFlat = 0.02;
HotPixThrSam = 0.02;
%% Body of function.
lambda    = EnergyConverter(EnergyPixelsizeDistance(1));
pixelsize = EnergyPixelsizeDistance(2);
distance  = EnergyPixelsizeDistance(3);
% Check input and output paths.
if InputPath(end) ~= '/'
    InputPath = [InputPath '/'];
end
outputPath = OutputPath;
if outputPath(end) ~= '/'
    outputPath = [outputPath '/'];
end
if ~exist(outputPath,'dir')
    mkdir(outputPath)
end
% Regularization parameters.
if size(alphaCTF_alphaTIE,2) == 1
    alphaCTF = alphaCTF_alphaTIE;
    alphaTIE = alphaCTF + log10(2*pi*lambda*distance/pixelsize^2);
else
    alphaCTF = alphaCTF_alphaTIE(1);
    alphaTIE = alphaCTF_alphaTIE(2);
end
fprintf('Regularization parameter: %.2f (CTF and TIE), %.2f (old TIE parameter)\n',alphaCTF,alphaTIE);
% Experimental parameters.
Energy    = EnergyPixelsizeDistance(1);
Pixelsize = EnergyPixelsizeDistance(2);
edpSize = numel(EnergyPixelsizeDistance);
switch edpSize
    case 3
        Dist1       = EnergyPixelsizeDistance(3);
        Dist2       = EnergyPixelsizeDistance(3);
    case 4
        Dist1       = EnergyPixelsizeDistance(3);
        Dist2       = EnergyPixelsizeDistance(4);
        fprintf('\nVarying distance measurement from %.3fm to %.3fm.\n',Dist1,Dist2)
end
%% Loop over different data sets (frog embryo stages).
tic
%% Create folders for retrieved-phase images will be stored.
% If folder already exists, MATLAB prints a warning.
if evalTIElo
    outputFolderTIElo = sprintf('TIElo_alpha%3.2f',alphaTIE);
    outputFolderTIElo = regexprep(outputFolderTIElo,'\.','p');
    if edpSize == 4
        outputFolderTIElo = [outputFolderTIElo '_distVar'];
    end
    outputPathTIElo   = [outputPath outputFolderTIElo '/'];
    if ~exist(outputPathTIElo,'dir')
        mkdir(outputPath,outputFolderTIElo);
    end
    
end
if evalTIEpnlo
    outputFolderTIEpnlo = sprintf('TIEpnlo_alpha%3.2f',alphaTIE);
    outputFolderTIEpnlo = regexprep(outputFolderTIEpnlo,'\.','p');
    if edpSize == 4
        outputFolderTIEpnlo = [outputFolderTIEpnlo '_distVar'];
    end
    outputPathTIEpnlo   = [outputPath outputFolderTIEpnlo '/'];
    if ~exist(outputPathTIEpnlo,'dir')
        mkdir(outputPath,outputFolderTIEpnlo);
    end
end
if evalCTF
    outputFolderCTF = sprintf('CTF_alpha%3.2f',alphaCTF);
    outputFolderCTF = regexprep(outputFolderCTF,'\.','p');
    if edpSize == 4
        outputFolderCTF = [outputFolderCTF '_distVar'];
    end
    outputPathCTF   = [outputPath outputFolderCTF '/'];
    if ~exist(outputPathCTF,'dir')
        mkdir(outputPath,outputFolderCTF);
    end
end
if BinaryFilterThreshold(1) > 0
    if isscalar(BinaryFilterThreshold)
        outputFolderCTFproj = sprintf('CTFproj_alpha%3.2f_binFilt%3.4f',alphaCTF,BinaryFilterThreshold);
    else
        outputFolderCTFproj = sprintf('CTFproj_alpha%3.2f_binFilt%3.4fblurW%02uS%0f',alphaCTF,BinaryFilterThreshold);
    end
    outputFolderCTFproj = regexprep(outputFolderCTFproj,'\.','p');
    if edpSize == 4
        outputFolderCTFproj = [outputFolderCTFproj '_distVar'];
    end
    outputPathCTFproj   = [outputPath outputFolderCTFproj '/'];
    if ~exist(outputPathCTFproj,'dir')
        mkdir(outputPath,outputFolderCTFproj);
    end
end
%% Loop over projections
loopcounter = 0;
% Determine data format.
DataFormatStruct = dir([InputPath FileNamePrefix '*']);
DataFormat = DataFormatStruct(1).name(end-2:end);
%% Read first image, set parameters and create the filter needed for phase retrieval.
intNameArray = dir([InputPath FileNamePrefix '*.*' DataFormat]);
FileName = [InputPath intNameArray(1).name];
switch lower(DataFormat)
    case {'edf'}
        gFT = pmedfread(FileName);
    case {'hdf'}
        gFT = hdfread(FileName, DataSetName, 'Index', IndexArray);
    case {'tif'}
        gFT = double(imread(FileName,DataFormat));
end
NumIm = numel(intNameArray);
%RefCell      = cell(1,60);
fprintf('Reading dark field images.')
for ff = 60:-1:1
    %fprintf('%2u ',ff)
    refindex = (NumIm-2*60)+ff;
    FileName = [InputPath intNameArray(refindex).name];
    im = double(imread(FileName,DataFormat));
    if ycrop > 0
        [dim1o dim2o] = size(gFT);
        ycroprange    = floor(dim2o/2)-ycrop/2+(1:ycrop);
        xcroprange    = floor(dim1o/2)-xcrop/2+(1:xcrop);
        im           = im(xcroprange,ycroprange);
    end
    im = FilterHotPixel(im,HotPixThrDark);
    RefCell{ff} = im;
end
fprintf('\n')
dark = median(cat(3,RefCell{:}),3);
domain(dark);
fprintf('Reading white field images.')
for ff = 60:-1:1
    %fprintf('%2u ',ff)
    refindex = (NumIm-1*60)+ff;
    FileName = [InputPath intNameArray(refindex).name];
    im = double(imread(FileName,DataFormat));
    if ycrop > 0
        %[dim1o dim2o] = size(gFT);
        ycroprange    = floor(dim2o/2)-ycrop/2+(1:ycrop);
        xcroprange    = floor(dim1o/2)-xcrop/2+(1:xcrop);
        im           = im(xcroprange,ycroprange);
    end
    im = FilterHotPixel(im-dark,HotPixThrFlat);
    RefCell{ff} = im;
end
fprintf('\n')
flat = median(cat(3,RefCell{:}),3);
%flat = flat - dark;
domain(flat);
PaddingFactor = Padding_FactorAndValue{1};
PaddingValue  = Padding_FactorAndValue{2};
% Wave length.
lambda        = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
% Crop image.
if ycrop > 0
    %[dim1o dim2o] = size(gFT);
    ycroprange    = floor(dim2o/2)-ycrop/2+(1:ycrop);
    xcroprange    = floor(dim1o/2)-xcrop/2+(1:xcrop);
    gFT           = gFT(xcroprange,ycroprange);
    fprintf('Cropping image!\nOriginal dimensions: %4i x %4i\n',dim1o,dim2o)
end
[dim1,dim2]   = size(gFT);
% Get new dimension whichs are the next-power-of-2 times PaddingFactor. The
% next-power-or-2 is done not to spoil computational power/precision because
% the fft of Matlab always pads to the next-power-of-2.
dimx          = PaddingFactor*2^nextpow2(dim1);
dimy          = PaddingFactor*2^nextpow2(dim2);
fprintf('Input dimensions: %4i x %4i\nProcessing dimensions: %4i x %4i\n',dim1,dim2,dimx,dimy)
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
% Number of projections.
%NumOfProj = numel(intNameArray);
NumOfProj = NumProj;
% Preallocate sinograms.
if doSino == 1
    SinoSliceNum = dim1/2;
elseif doSino > 1
    SinoSliceNum = doSino;
end
if doSino
    if evalTIElo
        sinotie = zeros(dim2,NumOfProj-ProjDecr);
    end
    if BinaryFilterThreshold(1) > 0
        sinopctf = zeros(dim2,NumOfProj-ProjDecr);
    end
end
%% Start looping.
fprintf('Projection which is currently processed and saved:\n')
for nn = StartProj:1:EndProj
    %% CTF filter are distance dependent!!
    % Filter for CTF phase retrieval.
    Distance    = Dist1 + (Dist2-Dist1)/(NumOfProj-1)*(nn-1);
    prefactor   = Pixelsize^2/(2*pi*lambda*Distance);
    SineXiEta   = sin(1/prefactor*(xi.^2 + eta.^2)/2);
    InverseSine = 1./(2*sign(SineXiEta).*(abs(SineXiEta))+10^-alphaCTF);
    % Projection filter for projected CTF phase retrieval.
    if BinaryFilterThreshold(1) > 0,
        % Create binary filter with threshold 'BinaryFilterThreshold(1)' for (1/prefactor*(xi.^2 + eta.^2)/2>pi/2)
        BinaryFilter = ones(dimx,dimy);
        BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
        if ~isscalar(BinaryFilterThreshold)
            hsize = BinaryFilterThreshold(2);
            sigma = BinaryFilterThreshold(3);
            BinaryFilter = imfilter(BinaryFilter,fspecial('gaussian',[hsize hsize],sigma));
            BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)/10) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
        end
    end
    %% Reading, padding and FT of data.
    FileName = [InputPath intNameArray(nn).name];
    switch lower(DataFormat)
        case {'edf'}
            gFT = pmedfread(FileName);
        case {'hdf'}
            gFT = hdfread(FileName, DataSetName, 'Index', IndexArray);
        case {'tif'}
            gFT = double(imread(FileName,DataFormat));
    end
    % Cropping.
    if ycrop > 0
            gFT = gFT(xcroprange,ycroprange);
    end 
    gFT = FilterHotPixel(gFT-dark,HotPixThrSam)./flat;
    % Pad intensity to dimensions [dimx dimy].
    gFT = padarray(gFT,[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PaddingValue,'pre');
    gFT = padarray(gFT,[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PaddingValue,'post');
    % Define intensity contrast g=I(x,y,z)/I(x,y,0)-1 and compute the Fourier
    % transform of g. For pure phase objects. <g> states the conservation of
    % flux <I>=1.  Substract mean is equivalent to the above definition of g in
    % the case of pure pahse objects. To play safe the mean is substracted
    % instead of 1. Since the mean value is not retrievable this doesn't change
    % the phase retrieval. The mean value (or (1,1) component of gFT) could
    % be such large that it's beyond numerical precision.
    gFT = gFT-mean(gFT(:));
    gFT = fft2(gFT);
    %% CTF PHASE RETRIEVAL
    % Divide FT[ctf] by the sine in Fourier space, invert
    % the Fourier transform, take real part (just to play safe, substract
    % mean (since it cannot be retrieved).
    if evalCTF
        ctf = real(ifft2(InverseSine.*gFT));
        ctf = ctf-mean(ctf(:));
        ctf = ctf(xcut,ycut);
    end
    %% Projected CTF
    if BinaryFilterThreshold(1) > 0
        % Projected CTF: Before inverse FT apply binary filter to FT of CTF
        % phase.
        ctfProjected  = real(ifft2(BinaryFilter.*InverseSine.*gFT));
        ctfProjected  = ctfProjected-mean(ctfProjected(:));
        ctfProjected  = ctfProjected(xcut,ycut);
        if doSino && (nn < NumOfProj -ProjDecr)
            sinopctf(:,nn) = ctfProjected(floor(SinoSliceNum),:);
        end
    end
    %% TIE PNLO (Perturbatively evaluated Next-to-Leading Order)
    % Fourier transform of TIE retrieved phase. Here, denotetd tieLO to safe
    % memory.
    if evalTIElo || evalTIEpnlo
        tieLO   = 1./(xi.^2 + eta.^2 + 10^-alphaTIE).*gFT;
    end
    % PNLO: Perturbatively evaluated Next-to-Leading Order correction.
    if evalTIEpnlo > 0
        % Compute the correction to the laplacian of the phase, denoted as
        % tieLOPNLO to save variable.
        phi_dx2 = ifft2(eta.^2.*tieLO);
        phi_dy2 = ifft2( xi.^2.*tieLO);
        tiePNLO = -real( ...
            + ifft2(eta.*tieLO).*( ifft2( (    (eta.^2 + xi.^2).*eta).*tieLO)) ...
            + ifft2( xi.*tieLO).*( ifft2( (xi.*(eta.^2 + xi.^2)     ).*tieLO)) ...
            + phi_dx2.^2 ...
            + phi_dx2.*phi_dy2 ...
            + phi_dy2.^2 ...
            + ifft2(xi.*eta.*tieLO).^2);
        % Substract the mean of the Laplacian of the NLO.
        tiePNLO = tiePNLO - mean(tiePNLO(:));
        % Compute the PNLO of the phase by means of Inversion of the
        % Laplacian using inverse Fourier transforms.
        tiePNLO = ifft2(1./(xi.^2 + eta.^2 + 10^-alphaTIE).*fft2(tiePNLO));
        % Take real part (just to play safe), rescale by prefactor, substract mean and crop back to
        % input dimensions.
        tiePNLO = prefactor*real(tiePNLO);
        tiePNLO = tiePNLO-mean(tiePNLO(:));
        % The PNLO phase map.
        tiePNLO = tiePNLO(xcut,ycut);
    end
    %% Linear TIE (Leading order).
    if evalTIElo
        % Compute TIE phase by means of Fourier
        % inversion. Take real part (to play safe).
        tieLO  = real(ifft2(tieLO));
        % Rescale by prefactor, substract mean and crop to
        % input dimensions.
        tieLO = tieLO - mean(tieLO(:));
        tieLO = prefactor*tieLO;
        tieLO = tieLO - mean(tieLO(:));
        % The LO-TIE phase map.
        tieLO = tieLO(xcut,ycut);
        if doSino && (nn < EndProj -ProjDecr)
            sinotie(:,1+nn-StartProj) = tieLO(floor(SinoSliceNum),:);
        end
    end
    if evalTIElo
        edfwrite(sprintf('%sphase_%04u.edf',outputPathTIElo',nn),tieLO','float32');
    end
    if evalTIEpnlo
        edfwrite(sprintf('%sphase_%04u.edf',outputPathTIEpnlo,nn),tiePNLO','float32');
    end
    if evalCTF
        edfwrite(sprintf('%sphase_%04u.edf',outputPathCTF,nn),ctf','float32');
    end
    if BinaryFilterThreshold(1) > 0
        edfwrite(sprintf('%sphase_%04u.edf',outputPathCTFproj,nn),ctfProjected','float32');
    end
    loopcounter  = loopcounter + 1;
    fprintf('%5u',nn)
    if mod(loopcounter,20) == 0
        fprintf('\n')
    end
end
if doSino
    if evalTIElo
        edfwrite(sprintf('%ssinoTIElin_slice%04u.edf',outputPath,SinoSliceNum),sinotie,'float32');
    end
    if BinaryFilterThreshold(1) > 0
        edfwrite(sprintf('%ssinoCTFproj_slice%04u.edf',outputPath,SinoSliceNum),sinopctf,'float32');
    end
end
fprintf('\n')
telapsed = toc;
fprintf('A total of %u images processed and saved in %gs (%gs/image).\n',loopcounter,telapsed,telapsed/loopcounter)
