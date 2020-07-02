function FlatCorOptXeno(doWrite,ParentPath,DataPrefix,FlatPrefix,DataSetsToProcess,RegionOfInterest,HotPixFiltThreshold,HotPixFiltPrintInfo,doReadRingCurrent)
% For details see function 'FlatCor'.

%% Default arguments.
if nargin < 1
    doWrite = 1;
end
if nargin < 2
    ParentPath = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo';
end
if nargin < 3
    DataPrefix = '*xeno';
end
if nargin < 4
    FlatPrefix = 'ref';
end
if nargin < 5
    DataSetsToProcess = 23;
end
if nargin < 6
    %first row: horizontal ROI, secod row: vertical ROI
      % mag 10x: 1.5 micron
      %RegionOfInterest{-1+1} = [ 1  600; 301 450]; sample out of field of
      %view
      RegionOfInterest{-1+2} = [201 1144; 1 746];
      RegionOfInterest{-1+3} = [201 1144; 1 746];
      RegionOfInterest{-1+4} = [221 1156; 67 780];
      RegionOfInterest{-1+5} = [221 1156; 65 794];
      RegionOfInterest{-1+6} = [208 1142; 49 728];
      RegionOfInterest{-1+7} = [208 1142; 49 778];
      RegionOfInterest{-1+8} = [202 1142; 21 728];
      RegionOfInterest{-1+9} = [204 1142; 1 674];
      RegionOfInterest{-1+10} = [192 1149; 20 742];
      % mag 4x: 3.75 micron
      RegionOfInterest{-1+11} = [97 500; 70 370];
      RegionOfInterest{-1+12} = [97 500; 87 404];
      RegionOfInterest{-1+13} = [91 500; 87 404];
      RegionOfInterest{-1+14} = [87 508; 171 470];
      % mag 10x: 1.5 micron, dist: 300 mm
      RegionOfInterest{-1+15} = [267 1308; 307 1068];
      RegionOfInterest{-1+16} = [267 1308; 307 1068];
      RegionOfInterest{-1+17} = [267 1308; 307 1068];
      RegionOfInterest{-1+18} = [267 1308; 307 1068];
      % mag 10x: 1.5 micron, dist: 500 mm
      RegionOfInterest{-1+19} = [247 1286; 287 1052];
      RegionOfInterest{-1+20} = [247 1286; 287 1052];
      RegionOfInterest{-1+21} = [247 1286; 287 1052];
      RegionOfInterest{-1+22} = [247 1286; 287 1052];
      % mag 20x: 0.75 micron, dist: 500 mm
      RegionOfInterest{-1+23} = [67 2006; 301 1768];
      RegionOfInterest{-1+24} = [67 2006; 301 1768];
      RegionOfInterest{-1+25} = [67 2006; 301 1768];
      % mag 10x: 1.5 micron, dist: 100 mm
      RegionOfInterest{-1+26} = [31 954; 15 750];
      %RegionOfInterest{27} = [31 954; 15 750];
      % mag 10x: 1.5 micron, dist: 400 mm
      RegionOfInterest{-1+27} = [1 961; 1 748];
      %RegionOfInterest{29} = [1 961; 1 748];
      % mag 20x: 0.75 micron, dist: 100 mm
      RegionOfInterest{-1+28} = [ 99  1962;345 1840];
      RegionOfInterest{-1+29} = [ 99  1962;345 1840];
      % mag 20x: 0.75 micron, dist: 400 mm
      RegionOfInterest{-1+30} = [ 78  1962;210 1678];
      RegionOfInterest{-1+31} = [ 78  1962;210 1678];
end
if nargin < 7
    HotPixFiltThreshold = 0.05;
end
if nargin < 8
    HotPixFiltPrintInfo = 0;
end
if nargin < 9
    doReadRingCurrent = 0;
end

%% Start looping over data sets and correcting data.
FlatCor(doWrite,ParentPath,DataPrefix,FlatPrefix,DataSetsToProcess,RegionOfInterest,HotPixFiltThreshold,HotPixFiltPrintInfo,doReadRingCurrent)

fprintf('\nFINISHED FLAT-/DARK-FIELD CORRECTION AND HOT-PIXEL-FILTERING OF THE DATA CONTAINED IN:\n %s\n',ParentPath);