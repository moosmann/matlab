% function header=writeheader(sizes,datatype,nbytes,nimage)
% used by edfwrite
% sizes is size(image to be written), cf. Dim_1 and Dim_2
% datatype is string, cf. DataType
% nbytes is number of bytes used for each written element
% nimage is probably number of the image, used for writing multiple images in a single file?
% origin: peter

function header=writeheader(sizes,datatype,nbytes,nimage)

headerlength=1024;

header=sprintf('%s\n','{');

header=[header headerstring('HeaderID','EH:000001:000000:000000','string')];
header=[header headerstring('Image',nimage,'integer')];
header=[header headerstring('ByteOrder','HighByteFirst','string')];
header=[header headerstring('DataType',datatype,'string')];
header=[header headerstring('Dim_1',sizes(1),'integer')];
header=[header headerstring('Dim_2',sizes(2),'integer')];
header=[header headerstring('Size',nbytes*sizes(1)*sizes(2),'integer')];
header=[header headerstring('Date',date,'string')];


header=[header blanks(headerlength-2-length(header)) sprintf('%s\n','}')];
