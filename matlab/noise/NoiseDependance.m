%% Noise analysis

for nn = 10:-1:1;
    [pha, int(:,:,nn) int0(:,:,nn)]=Lena(1.5,0,[20 2*nn/10 1e-6],4500);gn=int(:,:,nn)-int0(:,:,nn);
    fprintf('\nmean(int):%10g, 10^5mean(noise):%10g, mean(abs(noise)):%10g\n',mean(mean(int(:,:,nn))),mean(gn(:))*1e5,mean(abs(gn(:))));
    domain(int(:,:,nn))
end;
fprintf('\n')

disp(squeeze(max(max(int)))')
disp(squeeze(min(min(int)))')