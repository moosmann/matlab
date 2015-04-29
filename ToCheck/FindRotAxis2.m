function [RotAxis,mcor] = FindRotAxis2(ScanNum,DistNum,xmin,xmax)

if nargin<3,xmin=1;end;
if nargin<4,xmax=2048;end;


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
m1 = fliplr(flipud(rot90((edfread(sprintf('%s1500.edf',FilePrefix)) - dark)./f1500)));
m0 = m0(xmin:xmax,:);
m1 = m1(xmin:xmax,:);

mcor = real(ifft2(fft2(m1).*fft2(rot90(m0,2))));
%mcor = real(ifft(fft(m1,[],2).*fft(rot90(m0,2),[],2),[],2));

mcor = fftshift(mcor,2);

[maxpos,maxpos]=max(mcor,[],2);
RotAxis = (1024 - mean(maxpos))/2 + 1024;


fprintf(1,['RotAxis: ' ScanFolder ': ' num2str(RotAxis) '\n']);
