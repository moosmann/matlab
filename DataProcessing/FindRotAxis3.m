function [RotAxis,mcor] = FindRotAxis3(ScanNum,DistNum,xmin,xmax)

if nargin<3,xmin=1;end
if nargin<4,xmax=2048;end

ParentFolder = '/mnt/tomoraid3/tomo/ESRF_20100411_InhouseExperiment/';
ScanName     = {'CT_carbonMeshesPApet'; 'CT_graphiteMinePApet'; 'CT_meshesPApet'; 'CT_meshesPApet_B'};

ScanFolder   = [ParentFolder char(ScanName(ScanNum)) '/' char(ScanName(ScanNum)) ...
    '_' num2str(DistNum) '_/'];
FilePrefix   = [ScanFolder char(ScanName(ScanNum)) '_' num2str(DistNum) '_'];

% Read flat and dark fields.
dark  = edfread([ScanFolder 'darkend0000.edf']);
f0000 = edfread([ScanFolder 'refHST0000.edf']) - dark;
f1500 = edfread([ScanFolder 'refHST1500.edf']) - dark;

m0 = flipud(rot90((edfread(sprintf('%s0000.edf',FilePrefix)) - dark)./f0000));
m1 = rot90(rot90((edfread(sprintf('%s1500.edf',FilePrefix)) - dark)./f1500),2);
m0 = m0(xmin:xmax,:);
m1 = m1(xmin:xmax,:);

mcor = real(ifft2(fft2(m1).*fft2(rot90(m0,2))));

[~,maxpos]=max(mcor(:,1:1024),[],2);
meanmaxpos = mean(maxpos);
[~,maxpos2]=max(mcor(:,1025:end),[],2);
meanmaxpos2 = mean(maxpos2);

shift = (meanmaxpos + 1024 - meanmaxpos2)/4;
RotAxis = 1024 + shift;

fprintf(1,['RotAxis: ' ScanFolder ': ' num2str(RotAxis,'%8.8g') ' (' ...
    num2str(meanmaxpos) ' ' num2str(meanmaxpos2) ')\n']);

