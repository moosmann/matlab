% function count=edfwrite(filename,matrix,datatype)
% writes an image in esrf data format
% 
% datatype can be
%	uint8	= 8 bits	= UnsignedByte
%	uint16	= 16 bits	= UnsignedShort
%	uint32	= 32 bits	= UnsignedLong = UnsignedInteger
%	float32	= 		= Float = Real
%	float64 =		= DoubleValue = Double
% writes bigendian files

function edfwrite(filename,matrix,datatype,varargin)

switch nargin
    case 3
        head=[];
    case 4
        head=varargin{1};
end

switch datatype
	case 'uint8',
	esrfdatatype='UnsignedByte';
	nbytes=1;
	case 'uint16',
	esrfdatatype='UnsignedShort';
	nbytes=2;
	case 'uint32',
	esrfdatatype='UnsignedLong';
	nbytes=4;
	case 'float32',
	esrfdatatype='Float';
	nbytes=4;
	case 'float64',
	esrfdatatype='DoubleValue';
	nbytes=8;
end
nimage=1;

if isempty(head)
    head=writeedfheader(size(matrix),esrfdatatype,nbytes,nimage);
end

fid=fopen(filename,'w','b');
fwrite(fid,head);
fwrite(fid,matrix,datatype);
fclose(fid);
