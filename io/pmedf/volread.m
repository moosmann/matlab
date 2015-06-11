% function volume=volread(filename,format,skip,xsize,ysize,zsize,byteorder)
% reads a volume in raw format into the matrix image
%
% name of the file (a string)
% 
% other parameters are guessed if only the filename is given
% 
% format: datatype, see fread for possibilities
% skip: number of bytes skipped at start file
% xsize: first dimension (rows in Matlab)
% ysize: second dimension (columns in Matlab)
% zsize: third dimension (default is 1)
% byteorder: 'b' or 'l' (default); big endian or little endian 
% result is double
% origin: wolfgang and peter

function volume=volread(filename,format,skip,xsize,ysize,zsize,byteorder)

defaultformat = 'float32';
defaultzsize = 1;

computertype=computer;
if findstr(computertype,'linux')
	defaultbyteorder='l';
elseif findstr(computertype,'SOL')	
	defaultbyteorder='b';
else
	defaultbyteorder='l';
end

possibleskip=[0 210];


% if nargin < 5
% 	disp('At least 5 input arguments required:')
% 	help volread
% 	return
% end

guess = 0;

switch nargin
case 0
	help volread
	return
case 1
	guess = 1;
	byteorder=defaultbyteorder;
case 4
	ysize = xsize;
	zsize = defaultzsize;
	byteorder=defaultbyteorder;
case 5
	zsize = defaultzsize;
	byteorder=defaultbyteorder;
case 6
	byteorder=defaultbyteorder;
end  


fid=fopen(filename,'r',byteorder);

if fid == -1
	disp(sprintf('File %s does not exist!',filename))
else
	if guess
		switch defaultformat
		case 'float32'
			nbytes=4;
		end
		vfile=dir(filename);
		sfile=vfile.bytes;
		possize=sqrt((sfile-possibleskip)/nbytes);
		[a,b]=min(abs(possize-round(possize)));
		if (a==0)
			format=defaultformat;
			skip=possibleskip(b);
			xsize=possize(b);
			ysize=xsize;
			zsize=1;
			disp(sprintf('Size seems %d by %d',xsize,ysize))
		else
			zsize=0;
		end
	end % guess
	if zsize
		fseek(fid,skip,-1);
		volume=fread(fid,xsize*ysize*zsize,format);
		volume=reshape(volume,[xsize,ysize,zsize]);
	else
		disp('Problems estimating sizes, enter them manually.')
		help volread
	end
	fclose(fid);	
end
