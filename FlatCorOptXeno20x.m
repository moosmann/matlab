function FlatCorOptXeno20x(doWrite,ParentPath,DataPrefix,FlatPrefix,DataSetsToProcess,RegionOfInterest,HotPixFiltThreshold,HotPixFiltPrintInfo,doReadRingCurrent)
% For details see function 'FlatCor'.

%% Default arguments.
if nargin < 1
    doWrite = -400;
end
if nargin < 2
    ParentPath = '/mntdirect/_data_visitor/mi1057/bm05/pc/opt_xeno_stage12_tomo';
end
if nargin < 3
    DataPrefix = 'xeno';
end
if nargin < 4
    FlatPrefix = 'ref';
end
if nargin < 5
    DataSetsToProcess = [3 2 4];
end
if nargin < 6
    %first row: horizontal ROI, secod row: vertical ROI
    RegionOfInterest{1} = [ 99  1962;345 1840];
    RegionOfInterest{2} = RegionOfInterest{1};
    RegionOfInterest{3} = [ 105 1950; 201 1674];
    RegionOfInterest{4} = RegionOfInterest{3};
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
