function RecoLoopGPU(alpha,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,DataSetsToProcess,EnergyDistancePixelsize,ParentPath,DataFolderNamePrefix)
% See function 'RecoLoop' for detailed info.

%% Default parameters.
if nargin < 1
    alpha        = 2.5; %Regularization parameter
end
if nargin < 2
    evalTIElo    = 1;%Evaluate leading-order (lo) transport-of-intensity (TIE) expression TIElo
end
if nargin < 3
    evalTIEpnlo  = 0;%Evaluate perturbatively next-to-leading-order (nlo) TIE expression: TIEpnlo
end
if nargin < 4
    evalCTF      = 0;%Evaluate contrast-transfer-function (CTF) expression.
end
if nargin < 5
    BinaryFilterThreshold = [0.1 3 1];%if >0: evaluate projected CTF (PCTF), where value is the threshold for the projection filter.
end
if nargin < 6
    DataSetsToProcess = 0;
end
if nargin < 7
    EnergyDistancePixelsize = [20 0.94865 1.4e-6];%[Energy Distance Pixelsize] in metre.1
end
if nargin < 8
    ParentPath = '/mnt/tomoraid3/user/moosmann/Xenopus_ESRF_May2011/';
end
if nargin < 9
    DataFolderNamePrefix = 0;
end
Padding_FactorAndValue = {1 'symmetric'};

fprintf('\nReconstructing on GPUs\n');
%% Set/get parameters.
if ~iscell(EnergyDistancePixelsize)
    Energy    = EnergyDistancePixelsize(1);
    Distance  = EnergyDistancePixelsize(2);
    Pixelsize = EnergyDistancePixelsize(3);
end
%% Check path string ending.
if ParentPath(end) == '/'
    ParentPath = ParentPath(1:end-1);
end
%% Read folder names of all data sets.
if DataFolderNamePrefix == 0
    DataFolderNames      = dir(ParentPath);
    DataFolderNames(1:2) = [];
else
    DataFolderNames = dir([ParentPath '/data/' DataFolderNamePrefix '*']);
end
inputParentPath  = [ParentPath '/int/'];
outputPath = [ParentPath '/phase/'];
%% Loop over different data sets (frog embryo stages).
tic
totalcounter = 0;
if DataSetsToProcess == 0
    DataSetsToProcess = 1:numel(DataFolderNames);
end
for kk = DataSetsToProcess
    if iscell(EnergyDistancePixelsize)
        Energy    = EnergyDistancePixelsize{kk}(1);
        Distance  = EnergyDistancePixelsize{kk}(2);
        Pixelsize = EnergyDistancePixelsize{kk}(3);
    end
    % Create folders for phase images.
    if evalTIElo
        outputFolderTIElo = sprintf('%s/TIElo_alpha%3.2f',DataFolderNames(kk).name,alpha);
        outputFolderTIElo = regexprep(outputFolderTIElo,'\.','p');
        mkdir(outputPath,outputFolderTIElo);
        outputPathTIElo   = [outputPath outputFolderTIElo '/'];
    end
    if evalTIEpnlo
        outputFolderTIEpnlo = sprintf('%s/TIEpnlo_alpha%3.2f',DataFolderNames(kk).name,alpha);
        outputFolderTIEpnlo = regexprep(outputFolderTIEpnlo,'\.','p');
        mkdir(outputPath,outputFolderTIEpnlo);
        outputPathTIEpnlo   = [outputPath outputFolderTIEpnlo '/'];
    end
    if evalCTF
        outputFolderCTF = sprintf('%s/CTF_alpha%3.2f',DataFolderNames(kk).name,alpha);
        outputFolderCTF = regexprep(outputFolderCTF,'\.','p');
        mkdir(outputPath,outputFolderCTF);
        outputPathCTF   = [outputPath outputFolderCTF '/'];
    end
    if BinaryFilterThreshold(1) > 0
        if isscalar(BinaryFilterThreshold)
            outputFolderCTFproj = sprintf('%s/CTFproj_alpha%3.2f_binFilt%3.4f',DataFolderNames(kk).name,alpha,BinaryFilterThreshold);
        else
            outputFolderCTFproj = sprintf('%s/CTFproj_alpha%3.2f_binFilt%3.4fblurW%02uS%0f',DataFolderNames(kk).name,alpha,BinaryFilterThreshold);
        end
        outputFolderCTFproj = regexprep(outputFolderCTFproj,'\.','p');
        fprintf('\nOutput path: %s',outputPath);
        fprintf('\nOutput folder projected CTF: %s\n',outputFolderCTFproj);
        mkdir(outputPath,outputFolderCTFproj);
        outputPathCTFproj   = [outputPath outputFolderCTFproj '/'];
    end
    inputPath = [inputParentPath DataFolderNames(kk).name '/'];
    %% Loop over projections
    loopcounter = 0;
    %% Read first image, set paramters and create filters for phase retrieval.
    Intensity = pmedfread(sprintf('%sint_%04u.edf',inputPath,1));
    PaddingFactor = Padding_FactorAndValue{1};
    PaddingValue  = Padding_FactorAndValue{2};
    % Wave length.
    lambda        = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
    % Prefactor needed for TIE and CTF retrieval.
    prefactor     = Pixelsize^2/(2*pi*lambda*Distance);
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
        BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
        if ~isscalar(BinaryFilterThreshold)
            hsize = BinaryFilterThreshold(2);
            sigma = BinaryFilterThreshold(3);
            BinaryFilter = imfilter(BinaryFilter,fspecial('gaussian',[hsize hsize],sigma));
            BinaryFilter( (SineXiEta.^2<BinaryFilterThreshold(1)/10) & ((xi.^2 + eta.^2)>pi*prefactor) ) = 0;
        end
        BinaryFilter = gpuArray(BinaryFilter);
    end
    % Filter for linear TIE phase retrieval.
    InverseSine = 1./(2*sign(SineXiEta).*(abs(SineXiEta))+10^-alpha);
    InverseSine = gpuArray(InverseSine);
    xi          = gpuArray(xi);
    eta         = gpuArray(eta);
    fprintf('Projection of ''%s'' which is currently processed and saved:\n',DataFolderNames(kk).name)
    for nn =1:numel(dir([inputPath 'int_*.edf']))
        % Reading, padding and FT of data.
        Intensity = pmedfread(sprintf('%sint_%04u.edf',inputPath,nn));
        % Pad intensity to dimensions [dimx dimy].
        Intensity = padarray(Intensity,[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PaddingValue,'pre');
        Intensity = padarray(Intensity,[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PaddingValue,'post');
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
        %% CTF PHASE RETRIEVAL
        % Divide FT[ctf] by the sine in Fourier space, invert
        % the Fourier transform, take real part (just to play safe, substract
        % mean (since it cannot be retrieved).
        if evalCTF
            ctf       = gather(real(ifft2(InverseSine.*Intensity)));
            ctf       = ctf-mean(ctf(:));
            ctf       = ctf(xcut,ycut);
        end
        % Projected CTF
        if BinaryFilterThreshold(1) > 0
            % Projected CTF: Before inverse FT apply binary filter to FT of CTF
            % phase.
            ctfProjected  = gather(real(ifft2(BinaryFilter.*InverseSine.*Intensity)));
            ctfProjected  = ctfProjected-mean(ctfProjected(:));
            ctfProjected  = ctfProjected(xcut,ycut);
            % If BinaryFilter is needed for output uncomment next line.
            %BinaryFilter  = fftshift(gather(BinaryFilter));
        end
        %% TIE PNLO (Perturbatively evaluated Next-to-Leading Order)
        % Fourier transform of TIE retrieved phase. Here, denotetd tieLO to safe
        % memory.
        if evalTIElo || evalTIEpnlo
            tieLO       = 1./(xi.^2 + eta.^2 + 10^-alpha).*Intensity;
        end
        % PNLO: Perturbatively evaluated Next-to-Leading Order correction.
        if evalTIEpnlo > 0
            % Compute the correction to the laplacian of the phase, denoted as
            % tiePNLO to save variable.
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
            %tiePNLO = tiePNLO - mean(tiePNLO(:));
            % Compute the PNLO of the phase by means of Inversion of the
            % Laplacian using inverse Fourier transforms.
            tiePNLO = ifft2(1./(xi.^2 + eta.^2 + 10^-alpha).*fft2(tiePNLO));
            % Take real part (just to play safe), rescale by prefactor, substract mean and crop back to
            % input dimensions.
            tiePNLO = gather(prefactor*real(tiePNLO));
            tiePNLO = tiePNLO-mean(tiePNLO(:));
            % The PNLO phase map.
            tiePNLO = tiePNLO(xcut,ycut);
        end
        %% Linear TIE (Leading order).
        if evalTIElo
            % Compute TIE phase by means of Fourier
            % inversion. Take real part (to play safe).
            tieLO      = gather(real(ifft2(tieLO)));
            % Rescale by prefactor, substract mean and crop to
            % input dimensions.
            tieLO = tieLO - mean(tieLO(:));
            tieLO = prefactor*tieLO;
            tieLO = tieLO - mean(tieLO(:));
            % The LO-TIE phase map.
            tieLO = tieLO(xcut,ycut);
        end
        %% Save retrieved phase.
        if evalTIElo
            edfwrite(sprintf('%sphase_%04u.edf',outputPathTIElo,nn),tieLO,'float32');
        end
        if evalTIEpnlo
            edfwrite(sprintf('%sphase_%04u.edf',outputPathTIEpnlo,nn),tiePNLO,'float32');
        end
        if evalCTF
            edfwrite(sprintf('%sphase_%04u.edf',outputPathCTF,nn),ctf,'float32');
        end
        if BinaryFilterThreshold(1) > 0
            edfwrite(sprintf('%sphase_%04u.edf',outputPathCTFproj,nn),ctfProjected,'float32');
        end
        loopcounter  = loopcounter + 1;
        totalcounter = totalcounter + 1;
        fprintf('%5u',nn)
        if mod(loopcounter,20) == 0
            fprintf('\n')
        end
    end
    fprintf('\n')
end
telapsed = toc;
fprintf('A total of %u images processed and saved in %gs (%gs/image)\n',totalcounter,telapsed,telapsed/totalcounter)
