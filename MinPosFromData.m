function ftsum = MinPosFromData(DataSet,padfac,imcrement)
% Script to determine parameters from experimental data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    DataSet = 1;
end
if nargin < 2
    % Padding factor
    padfac = 1;
end
if nargin < 3
    imcrement = 1;
end
if nargin < 4
    SaveFile = '/mnt/tomoraid-LSDF/users/moosmann/matlabSave/ftsum.m';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Path to ESFR data
switch DataSet
    case 1
        dataPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell/int/Xenopus_4cell_20keV/';
    case 2
        dataPath = '/mnt/tomoraid-LSDF/tomo/ESRF_May2011_Xenopus/int/Xenopus_37Stage_Agar/';
end
% Get file names
intName = FilenameCell([dataPath 'int']);
% Image size for initialization and preallocation
[dim1,dim2] = size(pmedfread([dataPath intName{1}])');
% Initialization
ftsum = zeros(padfac*[dim1 dim2]);
% Loop over projections
imcount = 0;
for nn = 1:imcrement:numel(intName)
    fprintf(' %4u',nn);if mod(nn,30)==0;fprintf('\n');end
    % Read intensities
    im = pmedfread([dataPath intName{nn}])';
    %substract mean
    %im = im - mean(im(:));
    % FT of int and summation
    ftsum = ftsum + fftshift(fft2(padarray(im,(padfac-1).*[dim1 dim2]/2,'symmetric','both')));
    imcount = imcount + 1;
end
fprintf('\n')
ftsum = ftsum/imcount;
% Reduce images size of padded images
if padfac > 1
    ftsum = ftsum/length(ftsum(:));
    ftsum = ftsum(1:2:end,1:2:end)+ftsum(2:2:end,1:2:end)+ftsum(1:2:end,2:2:end)+ftsum(2:2:end,2:2:end);
    ftsum = ftsum*length(ftsum(:));
end
save(SaveFile,'ftsum')
domain(log(1+abs(ftsum)))
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

