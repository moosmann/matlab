function im=edfread(filename,varargin)
% edfread
% function im=edfread(filename,varargin)
%      reads an image in edf format into the matrix im
%      uses pmedf_read by Petr Mikulik if only filename is given
%      if both header and image are needed, use [hd,im]=pmedf_read(filename)
%      does not transpose the image, 1st dimension remains 1st dimension of original file
%      
%      arguments:
%      argument 1: filename
%      argument 2 (optional): vector of adjacent rows to be read
%      argument 3 (optional): vector of adjacent colums to be read
%      argument 4 (optional): vector of layers to be read
%      
%      examples:
%      im = edfread('dark.edf');
%          reads complete image
%      im = edfread('dark.edf',10:200,30:50);
%          reads subimage corresponding to row (Dim_1) 10 until 200 and column (Dim_2) 30 until 50
%      im = edfread('dark.edf',10:200,30:50,[1 4 5]);
%          reads layers 1, 4 and 5 of a 3-dimensional image (see Dim_3 and W. Ludwig)
%
%      See also: edfwrite, pmedf_read, pmedf_write

% Copyright (C) 2007 P. Cloetens
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

% Author: P. Cloetens <cloetens@esrf.fr>
%
% 2007-03-03 P. Cloetens <cloetens@esrf.fr>
% * Initial revision


usepm = 0;

switch nargin
    case 1
        usepm = exist('pmedf_read','file');
    case 3
        rows = varargin{1};
        columns = varargin{2};
    case 4
        rows = varargin{1};
        columns = varargin{2};
        layers =varargin{3};
end

if usepm
    [~,im] = pmedf_read(filename);
else
    fid=fopen(filename,'r');
    
    if fid==-1
        im=fid;
        fprintf(sprintf('!!! Error opening file %s !!!',filename))
    else
        hd=readedfheader(fid);
        
        byteorder=findheader(hd,'ByteOrder','string');
        fidpos=ftell(fid); % store present position in file
        fclose(fid);
        if strcmp(byteorder,'HighByteFirst')
            byteorder='b';
        else
            byteorder='l';
        end
        fid=fopen(filename,'r',byteorder); % re-open with good format
        fseek(fid,fidpos,0); % re-position at end of header
        
        xsize=findheader(hd,'Dim_1','integer');
        ysize=findheader(hd,'Dim_2','integer');
        zsize=findheader(hd,'Dim_3','integer');
        if isempty(zsize)
            zsize=1;
        end
        
        datatype=findheader(hd,'DataType','string');
        switch datatype
            case 'UnsignedByte',
                datatype='uint8';
                nbytes=1;
            case 'UnsignedShort',
                datatype='uint16';
                nbytes=2;
            case {'UnsignedInteger','UnsignedLong'}
                datatype='uint32';
                nbytes=4;
            case {'Float','FloatValue','FLOATVALUE','Real'}
                datatype='float32';
                nbytes=4;
            case 'DoubleValue'
                datatype='float64';
                nbytes=8;
                %etcetera
        end
        
        if isempty(who('rows'))
            rows=1:xsize;
        end
        if isempty(who('columns'))
            columns=1:ysize;
        end
        if isempty(who('layers'))
            layers=1:zsize;
        end
        
        if zsize==1
            fseek(fid,nbytes*(rows(1)-1+(columns(1)-1)*xsize),0);
            im=fread(fid,[length(rows),length(columns)],sprintf('%d*%s',length(rows),datatype),nbytes*(xsize-length(rows)));
        else
            j=1;
            for i=layers
                fseek(fid,1024+nbytes*(rows(1)-1+(columns(1)-1)*xsize+xsize*ysize*(i-1)),-1);
                im(:,:,j)=fread(fid,[length(rows),length(columns)],sprintf('%d*%s',length(rows),datatype),nbytes*(xsize-length(rows)));
                j=j+1;
            end
        end
        fclose(fid);
        
    end
end % if usepm
