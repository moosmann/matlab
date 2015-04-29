function NoiseToSignal(im)

%% Statistics
imMean   = mean(im,2);
imStd = std(im,0,2);
imSignal = imStd./imMean;
% Normlize signal
dim = size(imSignal,1);
x = round(dim/2) + (-round(dim/9):round(dim/9));
imSignalm = mean(imSignal(x));

%% Plot figure
figure
subplot(2,2,1),plot(imMean),title('Mean'),
subplot(2,2,2),plot(imStd),title('Standard Deviation'),
subplot(2,2,3),plot(imSignal),title('Noise to signal: StaDev/Mean')
subplot(2,2,4),plot(1:dim,imSignal/imSignalm,'blue',x(:),1,'red',1:dim,1.5,'green'),title('Normalized noise to signal: StaDev/Mean')