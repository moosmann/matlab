function filtered = CheckHotPixelFilterThresholds(dataFileNamePrefix,filterThresholds_DarkRefData)
% Read dark field image, first and last refHST image and first and last
% image of the sample. Apply hot-pixel filter to each image and compare
% printing.
% output: struct of filtered images

%% Default arguments.
if nargin < 1
    dataFileNamePrefix = 'Xenopus';
end
if nargin < 2
    filterThresholds_DarkRefData = [1.015 1.01 1.025];
end
%% Filter thresholds
filtThresDark   = filterThresholds_DarkRefData(1);
filtThresRef    = filterThresholds_DarkRefData(2);
filtThresSample = filterThresholds_DarkRefData(3);
%% Get name struct with corresponding file names.
dataFileNameStruct = dir([dataFileNamePrefix '*.edf']);
darkFileNameStruct = dir('dark*.edf');
refFileNameStruct  = dir('refHST*.edf');
%% Read dark and first and last reference images.
dark     = pmedfread(darkFileNameStruct.name)/15;
refFirst = pmedfread(refFileNameStruct(1).name);
refLast  = pmedfread(refFileNameStruct(end).name);
datFirst = pmedfread(dataFileNameStruct(1).name);
datLast  = pmedfread(dataFileNameStruct(end).name);
%% Apply hot-pixel filter.
% Dark field image.
fprintf(1,'DARK: %f\n',filtThresDark);
filtered.dark    = FilterHotPixel(dark,filtThresDark);
% Median filtered flat field images.
fprintf(1,'First refHST: %f\n',filtThresRef);
filtered.refFirst = FilterHotPixel(refFirst,filtThresRef);
fprintf(1,'Last refHST: %f\n',filtThresRef);
filtered.refLast = FilterHotPixel(refLast,filtThresRef);
% Sample projections.
fprintf(1,'First projection: %f\n',filtThresSample);
filtered.projFirst = FilterHotPixel(datFirst,filtThresSample);
fprintf(1,'Last projection: %f\n',filtThresSample);
filtered.projLast = FilterHotPixel(datLast,filtThresSample);
