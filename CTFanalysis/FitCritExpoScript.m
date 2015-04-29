%close all hidden
clear smax ii x CritExp CritExpFit
smax = 370:10:500;
smin = 360;
for ii = length(smax):-1:1
    x=smin:smax(ii);
    figure('Name',sprintf('Fit critical exponential: s_max = %f',smax(ii)))
    CritExpFit{ii} = FitCritExpo(x,cfMinPos(x));
    CritExp(ii) = CritExpFit{ii}.cexp;
end
figure('Name','Critical exponent VS s_smax')
plot(smax(:),CritExp(:),'-+')