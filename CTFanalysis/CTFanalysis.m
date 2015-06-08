%Analysis of the CTF model for increasing values of the absolute phase
%shift.
clear all

%% Parameters.
printPlots(1) = 1;
printPlots(2) = 0;
filestring = [ getenv('HOME'), '/data/test_pattern/lena/lena.tif' ];
blurring = [8 8 1];
energy = 10;
distance = 0.5;
pixelsize = 1.1e-6;
PhotonCountsForPoissonNoise = 0;
counts_min = PhotonCountsForPoissonNoise;
counts_max = PhotonCountsForPoissonNoise;
paddim = 1/1*1024;
SineArgPreFac = pi*EnergyConverter(energy)*distance/(paddim*pixelsize)^2;
rescalVec = 1:10:500;
MaxPhaseShift = 0.01*rescalVec;
numMaxPhaseShift = numel(rescalVec);
dt = 10/numMaxPhaseShift;
% Read image and pad it with zeros.
phase0 = zeros(paddim,paddim);
phase0(paddim/2+(-255:256),paddim/2+(-255:256)) = normat(double(imread(filestring)));
domain(phase0, 1, 'initial phase')
%% Blur image.
if blurring(1) > 0
    hsizex = blurring(1);
    hsizey = blurring(2);
    sigma  = blurring(3);
    phase0  = imfilter(phase0,fspecial('gaussian',[hsizex hsizey],sigma));
end;
domain(phase0, 1, 'blurred phase')
phase0 = normat(phase0);

tic
for kk = numMaxPhaseShift:-1:1
    phase = MaxPhaseShift(kk) * phase0;
    %% Compute intensity
    int = Propagation(phase,[energy distance pixelsize],1,'symmetric',0);
    if PhotonCountsForPoissonNoise > 0
        counts = PhotonCountsForPoissonNoise;
        if PhotonCountsForPoissonNoise*max(int(:))<2^16,
            int = imnoise(uint16(counts*int),'poisson');
            counts_min = min([int(:); counts_min]);
            counts_max = max([int(:); counts_max]);
            int = 1/counts*double(int);
        else
            fprintf(1,'Warning: Counts exceed 2^16\n');
        end;
    end;
    %domain(int, 1, 'Poisson-noised intensity')
    % Substract mean.
    int = int - mean(int(:));
    % Modulus of Fourier transform of intensity contrast
    fint = fftshift(fft2(int));
    afint = abs(fint);
    % Modulus of ratio of Fourier transform of intensity contrast to Fourier transform of object.
    %aratio = abs(fint./fftshift(fft2(phase)));
    %% Integration of absolut of FT of intensity over angular sector to crop truncation rod.
    % Loop over angular offsets.
    for offsetFac = 1
        % Define angular sectors.
        thetaoffset = offsetFac*0.5;
        thetastep = pi/paddim/2;
        theta = thetaoffset:thetastep:pi/2-thetaoffset;
        theta = cat(2,theta,theta+pi/2,theta+pi,theta+3*pi/2);
        % Loop along radial direction.
        for rr = paddim/2:-1:1
            numtheta = length(theta);
            xx = zeros(numtheta,1);
            %xxar = zeros(numtheta,1);
            for ii = 1:numtheta
                xx(ii) = afint(round(paddim/2+1+rr*cos(theta(ii))),round(paddim/2+1+rr*sin(theta(ii))));
                %xxar(ii) = aratio(round(paddim/2+rr*cos(theta(ii))),round(paddim/2+rr*sin(theta(ii))));
            end
            y(rr,kk) = sum(xx);%/length(theta);
            %yar(rr,kk) = sum(xxar)/length(theta);
        end
        % Plot integrated line.
        %figure(kk)
        %plot(1:length(y(:,kk)),y(:,kk))
    end
    %upperLimit =  MaxPhaseShift(kk)*1000;
    %lowerLimit = -MaxPhaseShift(kk)*100;
end

% sector mask
mm = zeros(size(phase0));
for rr = paddim/2:-1:1
    numtheta = length(theta);
    for ii = 1:numtheta
        mm(round(paddim/2+rr*cos(theta(ii))),round(paddim/2+rr*sin(theta(ii)))) = 1;
    end
end


% Normalise
for kk = 1:numMaxPhaseShift
    y(:, kk) = normat(y(:, kk));
end
fprintf('Time for computation of maps: %fs\n',toc)

%% Plot data
if 0
    figure('Name','Sequence of integrate Fourier space spectrum')
    xpfull = 1:size(y,1);
    for kk = 1:numMaxPhaseShift
        plot(xpfull, y(:, kk))
        legend(sprintf('S = %3u\n max phase = %5.3g', kk, MaxPhaseShift(kk)))
        pause(dt)
    end
end

xMin = 25;
x = xMin:180;

SinArg = SineArgPreFac*x.^2;
%% Fit integrated line.
tic
for kk = length(MaxPhaseShift):-1:1
    %figure('Name',sprintf('Maximal phase shift: %f',MaxPhaseShift(kk)))
    %cf{kk} = FitIntLine(x,y(x,kk),MaxPhaseShift(kk)*2,-MaxPhaseShift(kk)*1);
    %figure('Name',sprintf('Maximal phase shift (ratio): %f',MaxPhaseShift(kk)))
    %cfar{kk} = FitIntLine(x,yar(x,kk),MaxPhaseShift(kk)*1000,-MaxPhaseShift(kk)*100);
    %disp(cf{kk})
    %disp(cfar{kk})
    
    %cf{kk} = FitIntLineNoPlots(x,y(x,kk),MaxPhaseShift(kk)*2,-MaxPhaseShift(kk)*1);
    cf{kk} = FitIntLineNoPlots2(x(:), y(x, kk), SineArgPreFac);
    %cfar{kk} = FitIntLineNoPlots(x,yar(x,kk),MaxPhaseShift(kk)*1000,-MaxPhaseShift(kk)*100);
    
    % Determine values and positions aof minima and maxima.
    [cfMaxVal(kk), cfMaxPos(kk)] = max(cf{kk}(x));
    [cfMinVal(kk), cfMinPos(kk)] = min(cf{kk}(x));
    %[cfarMaxVal(kk), cfarMaxPos(kk)] = max(cfar{kk}(x));
    %[cfarMinVal(kk), cfarMinPos(kk)] = min(cfar{kk}(x));
end

fprintf('Time for computation of line fits: %fs\n',toc)
%% Values, positions, velocity and ratio of maxima and minima.
% Correct for cut (xMin) in radial direction.
cfMaxPos = cfMaxPos + xMin;
cfMinPos = cfMinPos + xMin;
% cfarMaxPos = cfarMaxPos + xMin;
% cfarMinPos = cfarMinPos + xMin;
% Compute velocity of the minima and maxima.
cfMaxValD1 = diff(cfMaxVal);
cfMaxPosD1 = diff(cfMaxPos);
cfMinValD1 = diff(cfMinVal);
cfMinPosD1 = diff(cfMinPos);
% cfarMaxValD1 = diff(cfarMaxVal);
% cfarMaxPosD1 = diff(cfarMaxPos);
% cfarMinValD1 = diff(cfarMinVal);
% cfarMinPosD1 = diff(cfarMinPos);
% Compute supporting points for velocity vectores.
rescalVecD1 = interp1(1:length(rescalVec),rescalVec,0.5+(1:(length(rescalVec)-1)));
% Substract rescaled ground-state (GS) value from Maximum
cfMaxValn   = cfMaxVal-cfMaxVal(1)*rescalVec;
cfMaxValnD1 = diff(cfMaxValn);
% cfarMaxValn   = cfarMaxVal-cfarMaxVal(1)*rescalVec;
% cfarMaxValnD1 = diff(cfarMaxValn);

%% Plot sequence
xpfull = 1:size(y,1);
figure('Name', 'Sequence of integrate Fourier space spectrum: data, data ROI, fit')
for kk = 1:numMaxPhaseShift    
    plot(xpfull, y(:, kk), '.g', ...
        x, y(x, kk), '.b', ...
        x, cf{kk}(x), '-r', ...
        cfMinPos(kk), cfMinVal(kk), 'om')
    legend(sprintf('S = %3u\n max phase = %5.3g', kk, MaxPhaseShift(kk)))
    pause(dt)
end

figure('Name','Minimum position VS rescaling')
plot(rescalVec,cfMinPos,'.')

    
%% Plot position and values of minimum and maxium.
if printPlots(2)
    figure('Name','Maximum position VS rescaling')
    plot(rescalVec,cfMaxPos,'.')
    figure('Name','Minimum position VS rescaling')
    plot(rescalVec,cfMinPos,'.')
    figure('Name','Maximum value VS rescaling')
    plot(rescalVec,cfMaxVal,'.')
    figure('Name','Minimum value VS rescaling')
    plot(rescalVec,cfMinVal,'.')
    figure('Name','Ratio of minimum to Maximum VS rescaling')
    plot(rescalVec,cfMinVal./cfMaxVal,'.')
    % Velocity plots.
    figure('Name','Velocity of maximum position VS rescaling')
    plot(rescalVecD1,cfMaxPosD1,'.')
    figure('Name','Velocity of maximum value VS rescaling')
    plot(rescalVecD1,cfMaxValD1,'.')
    figure('Name','Velocity of minimum position VS rescaling')
    plot(rescalVecD1,cfMinPosD1,'.')
    figure('Name','Velocity of minimum value VS rescaling')
    plot(rescalVecD1,cfMinValD1,'.')
    figure('Name','Ratio of velocities of Minima and Maxima VS rescaling')
    plot(rescalVecD1,cfMinValD1./cfMaxValD1,'.')
    % Plots of Rescaled GS subtracted maximum.
    figure('Name','Rescaled ground-state subtracted maximum value VS rescaling')
    plot(rescalVec,cfMaxValn,'.')
    figure('Name','Velocity of rescaled ground-state subtracted maximum value VS rescaling')
    plot(rescalVecD1,cfMaxValnD1,'.')
    figure('Name','(minimum)/(rescaled-GS-subtracted maximum) VS rescaling')
    plot(rescalVec,cfMinVal./cfMaxValn,'.')
    figure('Name','(velocitiy of minimum)/(velocity of rescaled-GS-subtracted maximum) VS rescaling')
    plot(rescalVecD1,cfMinValD1./cfMaxValnD1,'.')
    %% Plots for ratio
    % Plot position and values of minimum and maxium.
    figure('Name','Maximum position VS rescaling')
    plot(rescalVec,cfarMaxPos,'.')
    figure('Name','Minimum position VS rescaling')
    plot(rescalVec,cfarMinPos,'.')
    figure('Name','Maximum value VS rescaling')
    plot(rescalVec,cfarMaxVal,'.')
    figure('Name','Minimum value VS rescaling')
    plot(rescalVec,cfarMinVal,'.')
    figure('Name','Ratio of minimum to Maximum VS rescaling')
    plot(rescalVec,cfarMinVal./cfarMaxVal,'.')
    % Velocity plots.
    figure('Name','Velocity of maximum position VS rescaling')
    plot(rescalVecD1,cfarMaxPosD1,'.')
    figure('Name','Velocity of maximum value VS rescaling')
    plot(rescalVecD1,cfarMaxValD1,'.')
    figure('Name','Velocity of minimum position VS rescaling')
    plot(rescalVecD1,cfarMinPosD1,'.')
    figure('Name','Velocity of minimum value VS rescaling')
    plot(rescalVecD1,cfarMinValD1,'.')
    figure('Name','Ratio of velocities of Minima and Maxima VS rescaling')
    plot(rescalVecD1,cfarMinValD1./cfarMaxValD1,'.')
    % Plots of Rescaled GS subtracted maximum.
    figure('Name','Rescaled ground-state subtracted maximum value VS rescaling')
    plot(rescalVec,cfarMaxValn,'.')
    figure('Name','Velocity of rescaled ground-state subtracted maximum value VS rescaling')
    plot(rescalVecD1,cfarMaxValnD1,'.')
    figure('Name','(minimum)/(rescaled-GS-subtracted maximum) VS rescaling')
    plot(rescalVec,cfarMinVal./cfarMaxValn,'.')
    figure('Name','(velocitiy of minimum)/(velocity of rescaled-GS-subtracted maximum) VS rescaling')
    plot(rescalVecD1,cfarMinValD1./cfarMaxValnD1,'.')
end