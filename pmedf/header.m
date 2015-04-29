## Copyright (C) 2007 P. Cloetens
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

## header
## function ret = header( filename (, what (, kind)))
##      Read header or information within header of an EDF file
##      If second argument is given: reads value of 'what' from the headerstring
##
##      arguments:
##      argument 1: filename
##          name of te EDF-file
##      argument 2: what ( default: none, complete header )
##          possible values: any of the header lines or a motor/counter mnemonic  
##      argument 3: kind ( default: 'motor' )
##          possible values: 'integer', 'float', 'motor' or 'string'
##
##      examples:
##          totalheader = header('refHST0000.edf');
##          motorpos = header('darkend0000.edf','pmy');
##          Dim_1 = header('darkend0000.edf','Dim_1','integer');
##          pixel_size = header('darkend0000.edf','optic_used','float');
##          exposure_time = header('darkend0000.edf','count_time','float');
##
##      see also:
##          findheader, readedfheader, pmedf_read, pmedfread

## Author: P. Cloetens <cloetens@esrf.fr>
## 
## 2007-09-25 P. Cloetens <cloetens@esrf.fr>
## * Initial revision

function ret = header(filename, what, kind)

    switch nargin
        case 0
            help header
            return
        case 1
            totalheader = 1;
        case 2
            totalheader = 0;
            kind = 'motor';
        case 3
            totalheader = 0;
    endswitch

    fid = fopen(filename,'rb');
    if (fid != -1)
        headerstr = readedfheader(fid);
        if totalheader
	    ret = headerstr;
        else
	    ret = findheader(headerstr,what,kind);
        endif
        fclose(fid);
    else
        printf('The file %s can not be read\n',filename)
    endif

endfunction
