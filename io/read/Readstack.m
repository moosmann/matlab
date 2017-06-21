function stack = Readstack(InputPath,StepSize_or_VecOfImagesToRead,FilenamePattern)
% Read images (default: tif) into 3D stack.
%
% Written by Julian Moosmann, first version: 2010. long ago, last version:
% 2017-06-21

%% Default arguments.
if nargin < 1
    InputPath = '';
end
if nargin < 2
    StepSize_or_VecOfImagesToRead = 1;
end
if nargin < 3
    FilenamePattern = '*.tif';
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

%% Read filenames into cell
if isnumeric(InputPath)
    fprintf('\nERROR: First input argument must be string.\n\n')
    return
end
if isempty(InputPath)
    InputPath = '.';
end
CheckTrailingSlash(InputPath);
files = FilenameCell(sprintf('%s%s',InputPath,FilenamePattern));
NumFiles = numel(files);

%% Read images into stack
if numel(StepSize_or_VecOfImagesToRead) == 1
    filesToRead = 1:StepSize_or_VecOfImagesToRead:NumFiles;
else
    filesToRead = StepSize_or_VecOfImagesToRead;
end
NumFilesToRead = numel(filesToRead);
switch lower(FilenamePattern(end-2:end))
    case 'edf'
        for nn = NumFilesToRead:-1:1
            stack(:,:,nn) = pmedfread(sprintf('%s%s',InputPath,files{filesToRead(nn)}))';
        end
    case {'img', 'dar', 'ref', 'sli', 'sln'}
        for nn = NumFilesToRead:-1:1
            stack(:,:,nn) = read_dat_jm(sprintf('%s%s',InputPath,files{filesToRead(nn)}))';
        end
    otherwise        
        for nn = NumFilesToRead:-1:1
            stack(:,:,nn) = imread(sprintf('%s%s',InputPath,files{filesToRead(nn)}));
        end
end

%% Print info
fprintf('Found %u files matching string pattern ''%s'' in directory ''%s''\n',NumFiles,FilenamePattern,InputPath);
fprintf('Read %u into stack of dimension [%u %u %u] in %g s = %.2g min \n',NumFilesToRead,size(stack),toc,toc/60);