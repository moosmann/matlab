% calibrate_piezo.m
clear all;

inputPath = '/data/visitor/ma1135/id19/GI/alignment/piezo_cal';
inputPrefix = 'radio';
piezo_name = 'gridy';

imageSize = [2048 2048];
roi = [500 1500; 500 1500];
stepArray = 133:1:212;
gratingPeriod = 2.4; % in [um]

fitVals = zeros(length(stepArray), 2); % 1st: motor pos, 2nd: int val

for step = 1:length(stepArray)
    [header im] = pmedf_read(sprintf('%s/%s%4.4i.edf',inputPath,inputPrefix,stepArray(step)));
    if(step==1)
        imtool(im(roi(1,1):roi(1,2),roi(2,1):roi(2,2)));
    end
    motor.mne = pmedf_findInHeader(header, 'motor_mne', 'string');
    motor.pos = pmedf_findInHeader(header, 'motor_pos', 'string');
    piezo.pos = pmedf_findMnePos(motor.mne, piezo_name);
    piezo.value = pmedf_findPos(motor.pos, piezo.pos);
    fitVals(step,1)=piezo.value;
    fitVals(step,2)=mean(mean(im(roi(1,1):roi(1,2),roi(2,1):roi(2,2))));
end

plot(fitVals(:,1),fitVals(:,2));

dofit = true;
while(dofit)
    amp(1) = (max(fitVals(:,2),[],1)-min(fitVals(:,2),[],1))/2
    ampTemp = input('limits for amplitude estimation: ');
    amp(2) = amp(1)-ampTemp;
    amp(3) = amp(1)+ampTemp;
    ampoff(1) = (max(fitVals(:,2),[],1)+min(fitVals(:,2),[],1))/2;
    ampoff(2) = min(fitVals(:,2),[],1);
    ampoff(3) = max(fitVals(:,2),[],1)
    phaseoff(2) = -100;
    phaseoff(3) = 100;
    phaseoff(1) = input('estimate phase offset: ')
    period(1) = 1/input('estimate period in voltage interval: ');
    periodTemp = input('limits for period estimation: ');
    period(2) = 1/(1/period(1)+periodTemp);
    period(3) = 1/(1/period(1)-periodTemp)
    fitRes = fit_piezo(fitVals(:,1),fitVals(:,2),amp, ampoff,phaseoff,period);
    piezoPeriod = 1/fitRes.b;
    piezoSig = confint(fitRes,0.637);
    piezoPeriodSig = 1/fitRes.b^2*(piezoSig(2,2)-fitRes.b);
    fprintf('piezo voltage period: %.14f +/- %.14f\n ',piezoPeriod,piezoPeriodSig);
    fprintf('piezo volta per micron: %.14f +/- %.14f\n',piezoPeriod/gratingPeriod,gratingPeriod*piezoPeriodSig);
    dofit = input('fit again ? [1] [0]: ');
end
        
