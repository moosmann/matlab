## Copyright (C) 2006 M. Cloetens
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

## readvolinfo
##      reads a volinfo file (e.g. raw.info or vol.info) into an Octave structure
##      quick and dirty implementation
##      numerical data and strings should be correctly separated

## Author: P. Cloetens <cloetens@coral1>
## 
## 2006-09-01 P. Cloetens <cloetens@coral1>
## * Initial revision

function [ volinfostruct ] = readvolinfo (infofilename)
    if exist(infofilename,'file')
        fid = fopen(infofilename);
        
        while !isequal(tmp=fgetl(fid),-1)
            #disp(tmp)
            if (tmp(1) != '!') # this is not a command line
                eqindex = index(tmp,'=');
                if eqindex
                    eval(['volinfostruct.' deblank(tmp(1:eqindex-1)) '=' tmp(eqindex+1:end) ';'],['volinfostruct.' deblank(tmp(1:eqindex-1)) "='" deblank(tmp(eqindex+1:end)) "';"]);
                endif
            endif
        endwhile
        
        fclose(fid);
    else
        printf('%s does not exist\n',infofilename);
        volinfostruct = [];
    endif
endfunction
