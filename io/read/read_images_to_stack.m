function stack = read_images_to_stack(InputPath,StepSize_or_VecOfImagesToRead,FilenamePattern, raw_im_shape, parloop, verbose)
% Read images (default: tif) into 3D stack.
%
% Written by Julian Moosmann, first version: 2010. long ago, last version:
% 2017-06-21

%% Default arguments.
if nargin < 1
    InputPath = '.';
end
if nargin < 2
    StepSize_or_VecOfImagesToRead = 1;
end
if nargin < 3
    FilenamePattern = '*.tif*';
end
if nargin < 4
    raw_im_shape = 'kit';
end
if nargin < 5
    parloop = 0;
end
if nargin < 6
    verbose = 1;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

%% Read filenames into cell
if isnumeric(InputPath)
    fprintf('\nERROR: First input argument must be string.\n\n')
    return
end
CheckTrailingSlash(InputPath);
files = FilenameCell(sprintf('%s%s',InputPath,FilenamePattern));
NumFiles = numel(files);
if verbose
    fprintf('Found %u files matching string pattern ''%s'' in directory ''%s''\n',NumFiles,FilenamePattern,InputPath);
end
if NumFiles == 0
    error( 'No files found matching the pattern: %s', FilenamePattern )
end

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
    case 'raw'
        if strcmpi( raw_im_shape, 'kit' )
            raw_im_shape = [5120 3840];
        elseif strcmpi( raw_im_shape, 'ehd' )
            raw_im_shape = [3056 3056];
        end
        for nn = NumFilesToRead:-1:1
            filename = sprintf('%s%s',InputPath,files{filesToRead(nn)})';
            stack(:,:,nn) = read_raw( filename(1:end-4), raw_im_shape, 'uint16' );
        end
    otherwise        
        for nn = NumFilesToRead:-1:1
            stack(:,:,nn) = imread(sprintf('%s%s',InputPath,files{filesToRead(nn)}));
        end
end

%% Print info
if verbose
    fprintf('Read %u images into volume [%u %u %u] in %g s = %.2g min \n',NumFilesToRead,size(stack),toc,toc/60);
end
