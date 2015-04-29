close all;
tic
dataPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell/int/Xenopus_4cell_20keV/';
intName = FilenameCell([dataPath 'int']);
[dim1,dim2] = size(pmedfread([dataPath intName{1}])');
% Initialization
padfac = 2;
ftsum = zeros(padfac*[dim1 dim2]);
for nn = 1:400:numel(intName)
    fprintf(' %4u',nn);if mod(nn,30)==0;fprintf('\n');end
    im = pmedfread([dataPath intName{nn}])';
    %substract mean
    im = im - mean(im(:));
    ftsum = ftsum + fftshift(fft2(padarray(im,(padfac-1).*[dim1 dim2]/2,'symmetric','both')));
    ftsum = ftsum/length(ftsum(:));
end
if padfac > 1
    ftsum = ftsum(1:2:end,1:2:end)+ftsum(2:2:end,1:2:end)+ftsum(1:2:end,2:2:end)+ftsum(2:2:end,2:2:end);
end
fprintf('\n')
toc
ftsum2 = ftsum;
%itool(log(1+abs(ftsum)))

