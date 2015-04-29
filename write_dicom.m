function out = write_dicom( vol , filenamePrefix_str, BitDepth, MultiframeSingleFile)

%% Default Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    filenamePrefix_str = 'test';
end
if nargin < 3
    BitDepth = 16;
end
if nargin < 4
    MultiframeSingleFile = 0;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dimX, dimY, dimZ] = size( vol );

% info struct
info.Filename =  [filenamePrefix_str '.dcm'];
%info.FileModDate =  '18-Dec-2000 11:06:43';
info.FileSize =  525436;
info.Format =  'DICOM';
info.FormatVersion =  3;
info.Width =  dimX;
info.Height =  dimY;
info.BitDepth = BitDepth;
info.ColorType =  'grayscale';
%info.SelectedFrames =  [];
%info.FileStruct =  [1x1 struct];
%info.StartOfPixelData =  1140;
%info.FileMetaInformationGroupLength =  192;
%info.FileMetaInformationVersion =  [2x1 uint8];
%info.MediaStorageSOPClassUID =  '1.2.840.10008.5.1.4.1.1.7';
info.PhotometricInterpretation = 'MONOCHROME2';
%info.NumberOfFrames = dimX;
info.PixelSpacing = [1 1];
info.SliceThickness = 1;
info.IOD = 'CT Image Storage';
info.MultiframeSingleFile = MultiframeSingleFile;

switch BitDepth
    case 8
        fun = @ (x) uint8( x );
    case 16
        fun = @ (x) uint16( x );
end
    

if MultiframeSingleFile == 0
    CheckAndMakePath(filenamePrefix_str);
    filenamePrefix_str = [filenamePrefix_str '/' filenamePrefix_str];
end
filename = [filenamePrefix_str '.dcm'];

out = dicomwrite( reshape( fun( 2^BitDepth * normat( vol ) ) , [dimX dimY 1 dimZ]) , filename, info );
