function WriteImage(filenamePrefix,im,imFormat)
% Write input image 'im' to file 'filenamePrefix'.'imFormat'. By default
% 'imFormat' writes 'tif' images using 'write32bitTIF'. 'edf' uses
% 'edfrwrite', other formats use 'imwritesc'.
% If 'imFormat' can be a cell of formats.
%
% filenamePrefix: full path, but WITHOUT extension
% im: image to write in format speciefied by 'imFormat'
% imFormat: string. available format: 'tif' (default), 'edf', or any other
% format 'imwrite' can handle with. Write function: tif: write32bitTIF
% (D. Haenschke), edf: edfwrite (ESRF), other: imwritesc (ASTRA toolbox)
%
% Writen by Julian Moosmann, 2013-10-28
%
% WriteImage(filenamePrefix,im,imFormat)


%% Return if format is empty
if nargin == 3    
    if isempty(imFormat)
        return;
    end
end

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    if filenamePrefix(end-3) == '.'
        imFormat = filenamePrefix(end-2:end);
    else
        imFormat = 'tif';
    end
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(imFormat)
    if filenamePrefix(end-2:end) == imFormat
        filenamePrefix(end-3:end) = [];
    end
    filename = sprintf('%s.%s',filenamePrefix,imFormat);
    switch lower(imFormat)
        case {'tif'}
            write32bitTIF(filename,im);
        case 'edf'
            edfwrite(filename,im','float32');
        otherwise
            imwritesc(im,filename);
    end
elseif iscell(imFormat)
    for nn = 1:numel(imFormat)
        if filenamePrefix(end-2:end) == imFormat{nn}
            filenamePrefix(end-3:end) = [];
        end
        filename = sprintf('%s.%s',filenamePrefix,imFormat{nn});
        switch lower(imFormat{nn})
            case {'tif'}
                write32bitTIF(filename,im);
            case 'edf'
                edfwrite(filename,im','float32');
            otherwise
                imwritesc(im,filename);
        end
    end
end  