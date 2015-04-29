%Plot the ratio of the velocity of the maximum to the velocity of the
%minimum versus the scaling.
%If data is center and scaled, cftool says: x is normalized by mean 500.5 and std 288.8

ca

rescalVecMean =  mean(rescalVec);
rescalVecStd  = std(rescalVec);
xNorm = (rescalVec-rescalVecMean)/rescalVecStd;

zta = 50;
rescalVecPP = [-(zta-1):0 rescalVec];
cfMinValPP = [zeros(1,zta) cfMinVal];
rescalVecPPMean =  mean(rescalVecPP);
rescalVecPPStd  = std(rescalVecPP);
xPPnorm = (rescalVecPP-rescalVecPPMean)/rescalVecPPStd;

figure('Name','Maximum VS scaling: data and fit')
cfMaxValFit = Poly9Fit(rescalVec,cfMaxVal);
figure('Name','Minimum VS scaling: data and fit')
cfMinValFit = Poly9Fit(rescalVecPP,cfMinValPP);

for ii = 10:-1:1
    cfMaxValFitCof(ii) = eval(['cfMaxValFit.p' int2str(ii)]);
    cfMinValFitCof(ii) = eval(['cfMinValFit.p' int2str(ii)]);
end

cfMaxValFitD1Cof = polyder(cfMaxValFitCof);
cfMinValFitD1Cof = polyder(cfMinValFitCof);

figure('Name','Velocity of maximum: and polynomial from derivative of fit')
plot(rescalVec,polyval(cfMaxValFitD1Cof,xNorm))
figure('Name','Velocity of minimum: polynomial from derivative of fit')
plot(rescalVecPP,polyval(cfMinValFitD1Cof,xPPnorm))
xxx = 1:length(rescalVec);
ratioMaxMin = polyval(cfMaxValFitD1Cof,xNorm(xxx))./polyval(cfMinValFitD1Cof,xPPnorm(zta+xxx));
xxx = 1:floor(4/10*length(rescalVec));
% Create figure
figure1 = figure('XVisual',...
    '0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98404194812 29.67743169791],...
    'Name','Ratio of fitted velocities');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
box(axes1,'on');
hold(axes1,'all');
set(gcf,'position',[0 0, 1600 900])
set(gcf,'position',[0 0, 800 600])
set(0,'defaulttextinterpreter','none')
% Create plot
plot(rescalVec(xxx),ratioMaxMin(xxx),'LineWidth',2);
% Create label
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
saveas(gcf,'/home/moosmann/ol/RatioMaxMin4.eps','epsc2')
saveas(gcf,'/home/moosmann/ol/Fig-3.eps','epsc2')