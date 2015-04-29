function FlatCorXenopus4cell(doWrite,ParentPath,DataPrefix,FlatPrefix,DataSetsToProcess,RegionOfInterest,HotPixFiltThreshold,HotPixFiltPrintInfo,doReadRingCurrent)
% For details see function 'FlatCor'.

%% Default arguments.
if nargin < 1
    doWrite = -400;
end
if nargin < 2
    %ParentPath = '/mnt/tomoraid3/user/moosmann/Xenopus_4cell';
    %ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_Xenopus_4cell';
    ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell';
end
if nargin < 3
    DataPrefix = 'Xenopus';
end
if nargin < 4
    FlatPrefix = 'ref';
end
if nargin < 5
    DataSetsToProcess = 1;
end
if nargin < 6
    %first row: horizontal ROI, secod row: vertical ROI
    RegionOfInterest{1} = [   87 2010 ;  271 1810];
    RegionOfInterest{2} = [   87 2000 ;  251 1810];
end
if nargin < 7
    HotPixFiltThreshold = 0.03;
end
if nargin < 8
    HotPixFiltPrintInfo = 0;
end
if nargin < 9
    doReadRingCurrent = 0;
end

%% Start looping over data sets and correcting data.
FlatCor(doWrite,ParentPath,DataPrefix,FlatPrefix,DataSetsToProcess,RegionOfInterest,HotPixFiltThreshold,HotPixFiltPrintInfo,doReadRingCurrent)


