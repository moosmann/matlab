function 
close all;
tic
% Path to ESFR data
dataPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell/int/Xenopus_4cell_20keV/';
%dataPath = '/mnt/tomoraid-LSDF/tomo/ESRF_May2011_Xenopus/int/Xenopus_37Stage_Agar/';
% Get file names
intName = FilenameCell([dataPath 'int']);
% Image size for initialization and preallocation
[dim1,dim2] = size(pmedfread([dataPath intName{1}])');
% Padding factor
padfac = 2;
% Initialization
ftsum = zeros(padfac*[dim1 dim2]);
% Loop over projections
for nn = 1:400:numel(intName)
    fprintf(' %4u',nn);if mod(nn,30)==0;fprintf('\n');end
    im = pmedfread([dataPath intName{nn}])';
    %substract mean
    im = im - mean(im(:));
    ftsum = ftsum + fftshift(fft2(padarray(im,(padfac-1).*[dim1 dim2]/2,'symmetric','both')));
    ftsum = ftsum/length(ftsum(:));
end
% Reduce images size of padded images
if padfac > 1
    ftsum = ftsum(1:2:end,1:2:end)+ftsum(2:2:end,1:2:end)+ftsum(1:2:end,2:2:end)+ftsum(2:2:end,2:2:end);
end
fprintf('\n')
toc
ftsum2 = ftsum;
%itool(log(1+abs(ftsum)))

