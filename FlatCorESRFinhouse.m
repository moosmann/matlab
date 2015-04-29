function FlatCorESRFinhouse(doWrite,ParentPath,DataPrefix,FlatPrefix,DataSetsToProcess,RegionOfInterest,HotPixFiltThreshold,HotPixFiltPrintInfo,doReadRingCurrent)
% Process data sets in '/mnt/tomoraid3/user/moosmann/Xenopus_4cell'

%% Default arguments.
if nargin < 1
    doWrite = 0;
end
if nargin < 2
    ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_20100411_InhouseExperiment';
   
end
if nargin < 3
    DataPrefix = 'CT';
end
if nargin < 4
    FlatPrefix = 'ref';
end
if nargin < 5
    % set 0 to process all data sets
    DataSetsToProcess = 0;
end
if nargin < 6
    % set 0 for whole image
    RegionOfInterest = 0;
%     %first row: horizontal ROI, secod row: vertical ROI
%     % 24h_20keV_
%     RegionOfInterest{1} = [   1 2048 ;  421 1610];
%     % 2m_A_
%     RegionOfInterest{2} = [   91 1864 ;  1 2020];
%     % 2m_B_
%     RegionOfInterest{3} = [   265 1820 ; 165 1878  ];
%     % 7d_12keV_A_ !!SCAN WAS ABORTED, AGAROSE BOILED!!
%     %RegionOfInterest{4} = [   421 1586;  99 2048];
%     % 7d_20keV_
%     RegionOfInterest{4} = [ 333   1732; 241 1988];
%     % 7d_20keV_B_
%     RegionOfInterest{5} = [ 521 1648; 1 2048  ];
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

FlatCor(doWrite,ParentPath,DataPrefix,FlatPrefix,DataSetsToProcess,RegionOfInterest,HotPixFiltThreshold,HotPixFiltPrintInfo,doReadRingCurrent)


