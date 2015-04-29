function RecoLoopLCI(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,TomoSetsToProcess,EnergyDistancePixelsize,ParentPath,DataSet,InputFolderPrefix,PreProcessingFolder,LCIFilterWidth_CircleCross,doRingFilter)
% Main script for phase retrieval based on the script 'RecoLoop', but
% altered for 2-BM@APS beam time data. Calling script is intended to be
% RecoLoopAPS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default parameters.
if nargin < 1
    % Scalar: regularization parameter needed to fix the singularity encountered
    % when inverting the laplacian (TIE) or sine function (CTF).
    alphaCTF_alphaTIE = 2.5;%[TIE CTF] 
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
    evalCTF      = 0;
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
    % Scalar (0,default), or 1x?-Vector Data sets to process. Default (0)
    % is processes all data sets found under 'IntParentPath', else
    % scalar or 1x?-Vector of the data set number(s) that should be processed.
    TomoSetsToProcess = 1;%input(sprintf('Tomo sets to process (row vector): '));
end
if nargin < 7
    % 1x3-Vector, or cell of 1x3-Vectors: 1x3-Vector = [Energy Distance
    % Pixelsize] in metre. Can be cell array to specify for each data set
    % its own parameters. E.g.:
    EnergyDistancePixelsize = [30 0.620 2.2e-6];
end
if nargin < 8
    % String: trailing seperator ('/') not needed. Parent path containing
    % the subfolder where the intensity images are found over which will be
    % looped over. By default only subfolder starting with tomo are used. 
   ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/';
end
if nargin < 9 
   %DataSet = 'wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';  
   DataSet = 'wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';  
end
if nargin < 10
    % 0 (default), or string: part of prefix of data folder names. If set 0
    % the prefix is assumed to be the same as that of 'IntParentPath'.
    InputFolderPrefix = 'tomo';
end
if nargin < 11
    %PreProcessingFolder = 'int_filtLineSectionMedWidthH063V001_noCropping_filtSino';
    PreProcessingFolder = 'int_filtLineSectionMedWidthH063V001_noCropping';
end
if nargin < 12
    LCIFilterWidth_CircleCross = [0 0];
end
if nargin < 13
    doRingFilter = 1;
end
if ParentPath(end) ~= '/'
    ParentPath = [ParentPath '/'];
end
IntParentPath = [ParentPath 'int/' DataSet '/' PreProcessingFolder];
if IntParentPath(end) == '/'
    IntParentPath = IntParentPath(1:end-1);
end
Padding_ValueMethod = {1 'symmetric'};
if LCIFilterWidth_CircleCross(1) > 0
    CircularFilterWidth = LCIFilterWidth_CircleCross(1);
    CrossFilterWidth    = LCIFilterWidth_CircleCross(2);
    fprintf('Fourier space filtering of intensity maps. Inverse Gaussian width: [Circle Cross] = [%g %g]\n',LCIFilterWidth_CircleCross)
end
LCIFilterWidth_CircleCross_1 = LCIFilterWidth_CircleCross(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DoProjectedCTF = BinaryFilterThreshold(1);
%% Check for varying experimental parameters.
if ~iscell(EnergyDistancePixelsize)
    Energy    = EnergyDistancePixelsize(1);
    Distance  = EnergyDistancePixelsize(2);
    Pixelsize = EnergyDistancePixelsize(3);
    lambda    = EnergyConverter(Energy);
    % Regularization parameters.
    if size(alphaCTF_alphaTIE,2) == 1
        alphaCTF = alphaCTF_alphaTIE;
        alphaTIE = alphaCTF + log10(2*pi*lambda*Distance/Pixelsize^2);
    else
        alphaCTF = alphaCTF_alphaTIE(1);
        alphaTIE = alphaCTF_alphaTIE(2);
    end
end
%% Read folder names of all data sets.
if InputFolderPrefix == 0
    InputFolderNames      = dir([IntParentPath '/tomo*']);
    InputFolderNames(1:2) = [];
else
    InputFolderNames = dir([IntParentPath '/' InputFolderPrefix '*']);
end
%[ParentPath DataSet] = fileparts(IntParentPath);
IntParentPath = [IntParentPath '/'];
ParentPath = fileparts(ParentPath);
OutputPath = [ParentPath '/phase/' DataSet '/' PreProcessingFolder];
OutputPath = [OutputPath '/'];
% Print info.
fprintf('\nPHASE RETRIEVAL\nDATA SET: ''%s''',DataSet)
fprintf('\nPATH TO TOMOGRAMS: %s',IntParentPath)
fprintf('\nTOMOGRAMS TO PROCESS:')
if TomoSetsToProcess == 0
    fprintf(' ALL\n')
else
    fprintf(' %s\n',mat2str(TomoSetsToProcess))
end
%% Loop over different data sets (frog embryo stages).
tic;
totalcounter = 0;
if TomoSetsToProcess == 0
    TomoSetsToProcess = 1:numel(InputFolderNames);
end
if LCIFilterWidth_CircleCross(1) > 0
    FolderPostfix = sprintf('_FilterCircleW%03uCrossW%03u_Padding%u',CircularFilterWidth,CrossFilterWidth,Padding_ValueMethod{1});
else
    FolderPostfix = '';
end
if doRingFilter
    FolderPostfix = ['_filtRing' FolderPostfix];
end
for kk = TomoSetsToProcess
    if iscell(EnergyDistancePixelsize)
        Energy    = EnergyDistancePixelsize{kk}(1);
        Distance  = EnergyDistancePixelsize{kk}(2);
        Pixelsize = EnergyDistancePixelsize{kk}(3);
        lambda    = EnergyConverter(Energy);
        % Regularization parameters.
        if size(alphaCTF_alphaTIE,2) == 1
            alphaCTF = alphaCTF_alphaTIE;
            alphaTIE = alphaCTF + log10(2*pi*lambda*Distance/Pixelsize^2);
        else
            alphaCTF = alphaCTF_alphaTIE(1);
            alphaTIE = alphaCTF_alphaTIE(2);
        end
    end
    %% Create folders for retrieved-phase images will be stored.
    % If folder already exists, MATLAB prints a warning.
    warning('off','MATLAB:MKDIR:DirectoryExists');
    if evalTIElo
        outputFolderTIElo = sprintf('TIElo_alpha%3.2f%s/%s',alphaTIE,FolderPostfix,InputFolderNames(kk).name);
        outputFolderTIElo = regexprep(outputFolderTIElo,'\.','p');
        mkdir(OutputPath,outputFolderTIElo);
        OutputPathTIElo   = [OutputPath outputFolderTIElo '/'];
    end
    if evalTIEpnlo
        outputFolderTIEpnlo = sprintf('TIEpnlo_alpha%3.2f%s/%s',alphaTIE,FolderPostfix,InputFolderNames(kk).name);
        outputFolderTIEpnlo = regexprep(outputFolderTIEpnlo,'\.','p');
        mkdir(OutputPath,outputFolderTIEpnlo);
        OutputPathTIEpnlo   = [OutputPath outputFolderTIEpnlo '/'];
    end
    if evalCTF
        outputFolderCTF = sprintf('CTF_alpha%3.2f%s/%s',alphaCTF,FolderPostfix,InputFolderNames(kk).name);
        outputFolderCTF = regexprep(outputFolderCTF,'\.','p');
        mkdir(OutputPath,outputFolderCTF);
        OutputPathCTF   = [OutputPath outputFolderCTF '/'];
    end
    if DoProjectedCTF > 0
        if isscalar(BinaryFilterThreshold)
            outputFolderCTFproj = sprintf('CTFproj_alpha%3.2f_binFilt%3.4f%s/%s',alphaCTF,BinaryFilterThreshold,FolderPostfixInputFolderNames(kk).name);
        else
            outputFolderCTFproj = sprintf('CTFproj_alpha%3.2f_binFilt%3.4fblurW%02uS%0f%s/%s',alphaCTF,BinaryFilterThreshold,FolderPostfix,InputFolderNames(kk).name);
        end
        outputFolderCTFproj = regexprep(outputFolderCTFproj,'\.','p');
        mkdir(OutputPath,outputFolderCTFproj);
        OutputPathCTFproj   = [OutputPath outputFolderCTFproj '/'];
    end
    InputPath = [IntParentPath InputFolderNames(kk).name '/'];
    %% Loop over projections
    %loopcounter = 0;
    %% Read first image, set parameters and create the filter needed for phase retrieval.
    g_ft = pmedfread(sprintf('%sint_%04u.edf',InputPath,1));
    PadVal = Padding_ValueMethod{1};
    PadMeth  = Padding_ValueMethod{2};
    % Wave length.
    lambda        = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
    % Prefactor needed for TIE and CTF retrieval.
    prefactor     = Pixelsize^2/(2*pi*lambda*Distance);
    % Dimensions of input data array.
    [dim1,dim2]   = size(g_ft);
    % Get new dimension whichs are the next-power-of-2 times PadVal. The
    % next-power-or-2 is done not to spoil computational power/precision because
    % the fft of Matlab always pads to the next-power-of-2.
    dimx          = PadVal*dim1;%PadVal*2^nextpow2(dim1);
    dimy          = PadVal*dim2;%PadVal*2^nextpow2(dim2);
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
    % Filter for CTF phase retrieval.
    SineXiEta  = sin(1/prefactor*(xi.^2 + eta.^2)/2);
    InverseSine   = 1./(2*sign(SineXiEta).*(abs(SineXiEta))+10^-alphaCTF);
    % Projection filter for projected CTF phase retrieval.
    if DoProjectedCTF > 0,
        % Create binary filter with threshold 'BinaryFilterThreshold(1)' for (1/prefactor*(xi.^2 + eta.^2)/2>pi/2)
        BinaryFilter = ones(dimx,dimy);
        BinaryFilter( (SineXiEta.^2<mod(BinaryFilterThreshold(1),1)) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
        if ~isscalar(BinaryFilterThreshold)
            hsize = BinaryFilterThreshold(2);
            sigma = BinaryFilterThreshold(3);
            BinaryFilter = imfilter(BinaryFilter,fspecial('gaussian',[hsize hsize],sigma));
            BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)/10) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
        end
        if DoProjectedCTF > 1
            BinaryFilter( ((xi.^2+eta.^2)/prefactor/2/pi > floor(BinaryFilterThreshold(1))) ) = 0;
        end
    end
    % Artifact filter
    % Circular center-symmetric high-pass filter in Fourier space
    if sum(LCIFilterWidth_CircleCross) > 0
        PadFacX       = dimx/dim1;
        PadFacY       = dimy/dim2;
        LCIfilter = ones(size(xi));
        if CircularFilterWidth > 0
            LCIfilter = LCIfilter.*(1-exp(-((dimy*xi/PadFacY/CircularFilterWidth).^2+(dimx*eta/PadFacX/CircularFilterWidth).^2)/2));
        end
        if CrossFilterWidth > 0
            LCIfilter = LCIfilter.*(1-exp(-(dimy*xi/PadFacY/CrossFilterWidth).^2/2)).*(1-exp(-(dimx*eta/PadFacX/CrossFilterWidth).^2/2));
        end
    end
    %% Start looping.
    fprintf('PROCESSING TOMO SET: %s\n',InputFolderNames(kk).name)
    fprintf('INPUT PATH:  %s\n',InputPath)
    fprintf('OUTPUT PATH: %s\n',OutputPath)
    fprintf('SIZE OF PROJECTION: (horizontal x vertical) = %u x %u (unpadded), %u x %u (padded)\n',dim1,dim2,dimx,dimy)
    fprintf('REGULARIZATION PARAMETER: -log10(RegPar) = %g (CTF), %g (TIE)\n',alphaCTF,alphaTIE)
    fprintf('NUMBER OF PROJECTIONS PROCESSED:\n')
    InputNameStruct = dir([InputPath 'int_*.edf']);
    parfor nn = 1:numel(InputNameStruct)
        %% Reading, padding and FT of data.
        g_ft = pmedfread([InputPath InputNameStruct(nn).name]);
        % Pad intensity to dimensions [dimx dimy].
        g_ft        = padarray(g_ft,[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PadMeth,'pre');
        g_ft        = padarray(g_ft,[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PadMeth,'post');
        % Define intensity contrast g=I(x,y,z)/I(x,y,0)-1 and compute the Fourier
        % transform of g. For pure phase objects. <g> states the conservation of
        % flux <I>=1.  Substract mean is equivalent to the above definition of g in
        % the case of pure pahse objects. To play safe the mean is substracted
        % instead of 1. Since the mean value is not retrievable this doesn't change
        % the phase retrieval. The mean value (or (1,1) component of g_ft) could
        % be such large that it's beyond numerical precision.
        g_ft       = g_ft-mean(g_ft(:));
        g_ft       = fft2(g_ft);
        if doRingFilter
            g_ft1 = median(g_ft(:,[1:3 end-1:end]),2);
            g_ft(:,1) = g_ft1;
        end
        %% Apply filter
        if LCIFilterWidth_CircleCross_1 > 0   
            g_ft = LCIfilter.*g_ft;
        end
        %% CTF PHASE RETRIEVAL
        % Divide FT[ctf] by the sine in Fourier space, invert
        % the Fourier transform, take real part (just to play safe, substract
        % mean (since it cannot be retrieved).
        if evalCTF
            ctf       = real(ifft2(InverseSine.*g_ft));
            ctf       = ctf-mean(ctf(:));
            ctf       = ctf(xcut,ycut);
        end
        %% Projected CTF
        if DoProjectedCTF > 0
            % Projected CTF: Before inverse FT apply binary filter to FT of CTF
            % phase.
            ctfProjected  = real(ifft2(BinaryFilter.*InverseSine.*g_ft));
            ctfProjected  = ctfProjected-mean(ctfProjected(:));
            ctfProjected  = ctfProjected(xcut,ycut);
        end
        %% TIE PNLO (Perturbatively evaluated Next-to-Leading Order)
        % Fourier transform of TIE retrieved phase. Here, denotetd tieLO to safe
        % memory.
        if evalTIElo || evalTIEpnlo
            tieLO       = 1./(xi.^2 + eta.^2 + 10^-alphaTIE).*g_ft;
        end
        % PNLO: Perturbatively evaluated Next-to-Leading Order correction.
        if evalTIEpnlo > 0
            % Compute the correction to the laplacian of the phase, denoted as
            % tieLOPNLO to save variable.
            phi_dx2     = ifft2(eta.^2.*tieLO);
            phi_dy2     = ifft2( xi.^2.*tieLO);
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
            tieLO      = real(ifft2(tieLO));
            % Rescale by prefactor, substract mean and crop to
            % input dimensions.
            tieLO = tieLO - mean(tieLO(:));
            tieLO = prefactor*tieLO;
            tieLO = tieLO - mean(tieLO(:));
            % The LO-TIE phase map.
            tieLO = tieLO(xcut,ycut);
        end
        if evalTIElo
            edfwrite(sprintf('%sphase_%04u.edf',OutputPathTIElo,nn),tieLO,'float32');
        end
        if evalTIEpnlo
            edfwrite(sprintf('%sphase_%04u.edf',OutputPathTIEpnlo,nn),tiePNLO,'float32');
        end
        if evalCTF
            edfwrite(sprintf('%sphase_%04u.edf',OutputPathCTF,nn),ctf,'float32');
        end
        if DoProjectedCTF > 0
            edfwrite(sprintf('%sphase_%04u.edf',OutputPathCTFproj,nn),ctfProjected,'float32');
        end
        totalcounter = totalcounter + 1;
%         loopcounter  = loopcounter + 1;
%         fprintf('%5u',nn)
%         if mod(loopcounter,30) == 0
%             fprintf('\n')
%         end
    end
    fprintf('\n')
end
telapsed = toc;
fprintf('A total of %u images processed and saved in %g min, %g h, or (%gs/image)\n',totalcounter,telapsed/60,telapsed/3600,telapsed/totalcounter)
