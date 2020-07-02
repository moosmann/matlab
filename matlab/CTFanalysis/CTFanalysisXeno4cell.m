%Analysis of the CTF model for increasing values of the absolute phase
%shift.

clear all
%close all hidden

%% Parameters.
printPlots = 1;
%filestring = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell/int/Xenopus_4cell_20keV/int_0437.edf';
filestring = '/export/scratch1/moosmann/ESRF_MI1079_ID19_July2011_inlineTomo/int/Xenopus_4cell_20keV/int_tif/int_0800.tif';
fprintf('Intensity map used:\n%s\n',filestring)
int0 = imread(filestring);
fprintf('Size of intensity map: [%u x %u]\n',size(int0));
domain(int0)
edp = [20 0.945 0.745e-6];
int0Mean = mean(int0(:));
phase0 = real(ifft2(PhaseFilter('ctf',size(int0),edp,2.5,0.1,'double').*fft2(int0 - int0Mean)));
%[~, out] = Reco(int0,[5 2.5],[energy distance pixelsize],0,0,1,0.01);
%phase0 = out.ctfProjected;

%% Read image and pad it with zeros.
paddim = 2048;
SineArgPreFac = pi*EnergyConverter(edp(1))*edp(2)/(paddim*edp(3))^2;
phimax = 2;
steps = 2;
rescalVec = 1:(phimax-1)/(steps-1):phimax;
MaxPhaseShift = 1*rescalVec;
domain(phase0)
%ishow(phase0)
%phase0 = RemoveLowFreq(phase0,70);
%domain(phase0)
%fprintf('Size of retrieved phase map: [%u x %u]\n',size(phase0));
tic
for kk = length(MaxPhaseShift):-1:1
    phase = MaxPhaseShift(kk)*phase0;
    %% Compute intensity
    %int = Propagation(phase,[energy distance pixelsize],2,'symmetric',0); 
    int = Propagation2(phase,0,edp,0,1); 
    fprintf('Phase multiplicator: %g\n',MaxPhaseShift(kk))
    %domain(int)
    % Substract mean.
    int = int-mean(int(:));
    %% Modulus of Fourier transform of intensity contrast
    fint = fftshift(fft2(int));
    afint = abs(fint);
    %ishow(log(afint))
    %% Modulus of ratio of Fourier transform of intensity contrast to Fourier transform of object.
    aratio = abs(fint./fftshift(fft2(phase)));
    %% Integration of absolut of FT of intensity over angular sector which avoid truncation rod.
    % Loop over angular offsets.
    for offsetFac = 5
        % Define angular sectors.
        thetaoffset = offsetFac*0.05;
        thetastep = pi/2*paddim/2;
        theta = thetaoffset:thetastep:pi/2-thetaoffset;
        theta = cat(1,theta,theta+pi/2,theta+pi,theta+3*pi/2);
        % Loop over points along radial direction.
        for rr = paddim/2:-1:1
            for ii = length(theta):-1:1
                xx(ii) = afint(round(paddim/2+rr*cos(theta(ii))),round(paddim/2+rr*sin(theta(ii))));
                xxar(ii) = aratio(round(paddim/2+rr*cos(theta(ii))),round(paddim/2+rr*sin(theta(ii))));
            end
            y(rr,kk) = sum(xx)/length(theta);
            yar(rr,kk) = sum(xxar)/length(theta);
        end
        % Plot integrated line.
        %figure('Name',sprintf('phase shift multiplier: %g',MaxPhaseShift(kk)))
        %plot(1:length(y(:,kk)),y(:,kk))
    end
    %upperLimit =  MaxPhaseShift(kk)*1000;
    %lowerLimit = -MaxPhaseShift(kk)*100; 
end
fprintf('Time for computation of maps: %fs\n',toc)
xMin = 150;
xMax = 240;
x = xMin:xMax;
SinArg = SineArgPreFac*x.^2;
%xMin = 170;
%x = xMin:230;
%% Fit integrated line.
tic
for kk = length(MaxPhaseShift):-1:1
    % Fit data points.
    cf{kk} = FitIntLineXeno4cell(double(x),double(y(x,kk)),SineArgPreFac,MaxPhaseShift(kk));
    % Determine values and positions aof minima and maxima.
    [cfMaxVal(kk), cfMaxPos(kk)] = max(cf{kk}(x));
    [cfMinVal(kk), cfMinPos(kk)] = min(cf{kk}(x));
end
fprintf('Time for computation of line fits: %fs\n',toc)
%% Values, positions, velocity and ratio of maxima and minima.
% Correct for cut (xMin) in radial direction.
cfMaxPos = cfMaxPos + xMin;
cfMinPos = cfMinPos + xMin;
fprintf('Maximum position: ')
fprintf('%u, ',cfMaxPos)
fprintf('\nMinimum position: ')
fprintf('%u, ',cfMinPos)
fprintf('\n')
% cfarMaxPos = cfarMaxPos + xMin;
% cfarMinPos = cfarMinPos + xMin;
% Compute velocity of the minima and maxima.
cfMaxValD1 = diff(cfMaxVal);
cfMaxPosD1 = diff(cfMaxPos);
cfMinValD1 = diff(cfMinVal);
cfMinPosD1 = diff(cfMinPos);

figure('Name','Data and fits vs index')
plot(x,cat(2,cf{1}(x),cf{2}(x),y(x,:)),'.')

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
% %% Plot position and values of minimum and maxium.
% if 0
%     figure('Name','Maximum position VS rescaling')
%     plot(rescalVec,cfMaxPos,'.')
%     figure('Name','Minimum position VS rescaling')
%     plot(rescalVec,cfMinPos,'.')
%     figure('Name','Maximum value VS rescaling')
%     plot(rescalVec,cfMaxVal,'.')
%     figure('Name','Minimum value VS rescaling')
%     plot(rescalVec,cfMinVal,'.')
%     figure('Name','Ratio of minimum to Maximum VS rescaling')
%     plot(rescalVec,cfMinVal./cfMaxVal,'.')
%     % Velocity plots.
%     figure('Name','Velocity of maximum position VS rescaling')
%     plot(rescalVecD1,cfMaxPosD1,'.')
%     figure('Name','Velocity of maximum value VS rescaling')
%     plot(rescalVecD1,cfMaxValD1,'.')
%     figure('Name','Velocity of minimum position VS rescaling')
%     plot(rescalVecD1,cfMinPosD1,'.')
%     figure('Name','Velocity of minimum value VS rescaling')
%     plot(rescalVecD1,cfMinValD1,'.')
%     figure('Name','Ratio of velocities of Minima and Maxima VS rescaling')
%     plot(rescalVecD1,cfMinValD1./cfMaxValD1,'.')
%     % Plots of Rescaled GS subtracted maximum.
%     figure('Name','Rescaled ground-state subtracted maximum value VS rescaling')
%     plot(rescalVec,cfMaxValn,'.')
%     figure('Name','Velocity of rescaled ground-state subtracted maximum value VS rescaling')
%     plot(rescalVecD1,cfMaxValnD1,'.')
%     figure('Name','(minimum)/(rescaled-GS-subtracted maximum) VS rescaling')
%     plot(rescalVec,cfMinVal./cfMaxValn,'.')
%     figure('Name','(velocitiy of minimum)/(velocity of rescaled-GS-subtracted maximum) VS rescaling')
%     plot(rescalVecD1,cfMinValD1./cfMaxValnD1,'.')
% %     %% Plots for ratio
% %     % Plot position and values of minimum and maxium.
% %     figure('Name','Maximum position VS rescaling')
% %     plot(rescalVec,cfarMaxPos,'.')
% %     figure('Name','Minimum position VS rescaling')
% %     plot(rescalVec,cfarMinPos,'.')
% %     figure('Name','Maximum value VS rescaling')
% %     plot(rescalVec,cfarMaxVal,'.')
% %     figure('Name','Minimum value VS rescaling')
% %     plot(rescalVec,cfarMinVal,'.')
% %     figure('Name','Ratio of minimum to Maximum VS rescaling')
% %     plot(rescalVec,cfarMinVal./cfarMaxVal,'.')
% %     % Velocity plots.
% %     figure('Name','Velocity of maximum position VS rescaling')
% %     plot(rescalVecD1,cfarMaxPosD1,'.')
% %     figure('Name','Velocity of maximum value VS rescaling')
% %     plot(rescalVecD1,cfarMaxValD1,'.')
% %     figure('Name','Velocity of minimum position VS rescaling')
% %     plot(rescalVecD1,cfarMinPosD1,'.')
% %     figure('Name','Velocity of minimum value VS rescaling')
% %     plot(rescalVecD1,cfarMinValD1,'.')
% %     figure('Name','Ratio of velocities of Minima and Maxima VS rescaling')
% %     plot(rescalVecD1,cfarMinValD1./cfarMaxValD1,'.')
% %     % Plots of Rescaled GS subtracted maximum.
% %     figure('Name','Rescaled ground-state subtracted maximum value VS rescaling')
% %     plot(rescalVec,cfarMaxValn,'.')
% %     figure('Name','Velocity of rescaled ground-state subtracted maximum value VS rescaling')
% %     plot(rescalVecD1,cfarMaxValnD1,'.')
% %     figure('Name','(minimum)/(rescaled-GS-subtracted maximum) VS rescaling')
% %     plot(rescalVec,cfarMinVal./cfarMaxValn,'.')
% %     figure('Name','(velocitiy of minimum)/(velocity of rescaled-GS-subtracted maximum) VS rescaling')
% %     plot(rescalVecD1,cfarMinValD1./cfarMaxValnD1,'.')
% end