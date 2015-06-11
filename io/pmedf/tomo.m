%% Copyright (C) 2006 P. Cloetens
%% 
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%% 
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

%% tomo
%% function []=tomo(varargin)
%%      used to create 	median/average reference files (refHST*.edf)
%%			average dark file (dark.edf)
%%      interactive if called without arguments
%%      
%%      arguments:
%%      argument 1: directory to be processed, path is absolute or relative ( default: current directory )
%%      argument 2: way images are 'summed': median or average ( default: 'm' )
%%      argument 3: value of save_ram ( default 2 )
%%                  save_ram = 0 : process complete images at once
%%                  save_ram = 1 : process images line by line
%%                  save_ram = 2 : call PyHST for processing
%%
%% origin: Wolfgang & Peter & Lukas
%% Date: 08/04/2003
%% Version: 0.1

%% Author: P. Cloetens <cloetens@esrf.fr>
%% 
%% 2006-07-28 P. Cloetens <cloetens@esrf.fr>
%% * Initial revision
%% * Changes with respect to archive version:
%% * creates dark.edf if not existing
%% * checks more correctly for correct file size of refHST files
%% * allow for relative paths
%% 2007-05-04 PC
%% * adapt to case where nvue is not a multiple of refon (360 degrees)
%% * last reference can be either nvue or refon*round(nvue/refon)
%% * possibility to call tomo with 1, 2 or 3 arguments
%% * possibility to calculate the median with PyHST (DO_PROJECTION_MEDIAN) 
%% * for dark.edf, save d+0.5 (edfwrite truncates in conversion to uint16)
%% * check if reference needs to be created for all values of save_ram
%% * display of saved reference is refreshed instead of multiple lines
%% 2007-05-10 PC
%% * make save_ram = 2 the default (use PyHST)
%% * save ref+0.5 as conversion to uint16 truncates
%% * remove parasitic PyHST file ref.vol
%% * remove global variable FT_DEVSERVER
%% * make refsname external function with 3 arguments
%% 2007-05-21 PC
%% * foresee arbitrary indices for reference list
%% * indices 0:refon:nvue and nvue always included

function []=tomo(varargin)

global FT_DELETESINGLEREF

if isempty(FT_DELETESINGLEREF)
    deletesingleref = 1;
else
    deletesingleref = FT_DELETESINGLEREF;
end


refroot='refHST';

default_summed = 'm';
default_save_ram = 2;
% save_ram = 0 : process complete images at once
% save_ram = 1 : process images line by line
% save_ram = 2 : call PyHST for processing
summed = [];
save_ram = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Information on job   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
	% interactive information
	% n1: directory, default is present directory
	n1=cleandirectoryname;
	% summed: option for making refHST-files, default is median
	summed = default_summed;
	summed = inputwdefault('(a)verage or (m)edian of reference images',summed);
	if summed=='m'
		disp('Median of reference images is used')
	else
		disp('Average of reference images is used')
	end
    case 1
	n1 = varargin{1};
	summed = [];
    case 2
	n1 = varargin{1};
	summed = varargin{2};
    case 3
	n1 = varargin{1};
	summed = varargin{2};
        save_ram = varargin{3};
endswitch

% just in case relative path is given
if isempty(n1)
    n1 = cleandirectoryname(pwd);
elseif n1(1) ~= '/'
    n1 = sprintf('%s/%s',cleandirectoryname(pwd),n1);
else
    n1 = cleandirectoryname(n1);
end
if isempty(summed)
    summed = default_summed;
endif
if isempty(save_ram)
    save_ram = default_save_ram;
endif

% n2: fileprefix, taken from directory name
pos=findstr(n1,'/');
n2=n1(pos(end)+1:end);

% information from edf-files

hd = header([n1 '/' n2 '0000.edf']);
X=findheader(hd,'Dim_1','integer');
Y=findheader(hd,'Dim_2','integer');

% information from info-file

fp=fopen([n1 '/' n2,'.info'],'r');
if fp~=-1 % *.info exists
    hd=fscanf(fp,'%c');
    fclose(fp);
    nvue =findheader(hd,'TOMO_N','integer');
    refon=findheader(hd,'REF_ON','integer');
    nref =findheader(hd,'REF_N','integer');
    ndark=findheader(hd,'DARK_N','integer');
else	% *.info does not exist
    printf('%s.info could not be found! Please enter parameters:\n',n2);
    n2   =inputwdefault('prefix',n2);
    nvue =input('number of projections        : ');
    refon=input('refon ?                      : ');
    nref =input('number of refs ?             : ');
    ndark=input('number of darks ?            : ');
end

% determine if device server was used
if ~isempty(refstobeconverted = glob('ref1????.edf'))
    % refs name convention is ref10000.edf etc, some files remain in place
    devserver = 0;
elseif ~isempty(refstobeconverted = glob('ref0000_????.edf'))
    % refs name convention is ref0000_0000.edf etc, some files remain in place 
    devserver = 1;
else
    % nothing to be converted
    devserver = 1;
        
end
% create vector with index of references to be converted
if ~isempty(refstobeconverted)
    index_refstobeconverted = zeros(1,length(refstobeconverted));
    for k = 1:length(index_refstobeconverted)
        pos = findstr(refstobeconverted{k},'.edf');
        index_refstobeconverted(k) = str2num(refstobeconverted{k}(pos-4:pos-1));
    end
    printf('Found %d series of references to be converted/removed\n',length(index_refstobeconverted));
else
    index_refstobeconverted = [];
endif
datatype = 'uint16';
if (~devserver) & (save_ram == 2)
    save_ram = 1;
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Creating dark file   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist([n1 '/dark.edf'],'file')
    disp('Creating dark.edf ...')
    if isempty(ndark)
        ndark=header([n1 'darkend0000.edf'],'no_summed','integer');
        if isempty(ndark)
            disp('Number of dark images is not specified, exiting !!!')
            return
        end
    end
    d=edfread([n1 '/darkend0000.edf'])/ndark;
    if isequal(size(d),[X Y])
        edfwrite([n1 '/dark.edf'],d+0.5,datatype);
    else
        disp('Dark does not exist or size is wrong')
        disp('Please wait until end of scan')
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Creating reference files   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we try to intercept the case where nvue is not a multiple of refon (360 degrees)
index_references = unique([0:refon:nvue nvue]);
% this is the usual list of indices for references
% maybe more references were taken, we concatenate with index_refstobeconverted
number_references = length(index_references);
index_references = unique([index_references index_refstobeconverted]);
if (length(index_references) > number_references)
    printf('More references than expected, all will be converted\n')
endif
sizenoheader = 2*X*Y;
firstconversion = 1;

if summed=='m' % use median
    switch save_ram
    case 1
    % process line by line
        f=zeros(X,nref);
        ref=zeros(X,Y);
        % determination byteorder of reference files
        if isempty(index_refstobeconverted)
            % no single reference files available, probably everything is done
            byteorder='l';
        else	
            byteorder=header([n1 '/' refsname(1,index_refstobeconverted(1),devserver)],'ByteOrder','string');
            if strcmp(byteorder,'HighByteFirst')
                    byteorder='b';
            else
                    byteorder='l';
            endif
        endif

        for n = index_references
            refname=sprintf('%s/%s%4.4i.edf',n1,refroot,n);

            % determine if reference needs to be created
            createref = ~checkreference(refname,sizenoheader);
            
            if createref  % reference does not exist or has wrong size and will be created
                if firstconversion
                    disp('Creating references ...')
                    firstconversion = 0;
                endif
	        % open the single reference files and move to end of header
                if isempty(find(index_refstobeconverted == n))
                    n_used = refon*round(n/refon);
                else
                    n_used = n;
                endif
                for k=1:nref 
		    name=[n1 '/' refsname(k,n_used,devserver)];
		    fid(k)=fopen(name,'r',byteorder);
		    hd=readedfheader(fid(k));
	        endfor
                
	        % read line by line and calculate median
	        for l = 1:Y
		    for k=1:nref
			    f(:,k)=fread(fid(k),X,datatype);
		    endfor
		    ref(:,l)=median(f,2);
	        endfor

	        % close the single reference files
	        for k = 1:nref
		    fclose(fid(k));
	        end

    	        % store the reference file
	        edfwrite(refname,ref+0.5,datatype);
   	        printf('Reference %s%4.4i.edf written\r',refroot,n)

                % in case last reference is 999 we also create a link 1000 to this file
                if (n == nvue)
                    nround = refon*round(n/refon);
                    refnameround = sprintf('%s/%s%4.4i.edf',n1,refroot,nround);
                    if (n ~= nround) & (~exist(refnameround,'file'))
                        link(refname,refnameround);
                    endif
                endif
            endif % creation new ref file
        endfor
    case 2
    % processing done by PyHST
        refpyhstname = 'ref_median.edf';
        refpyhstparfilename = 'ref.par';
        for n = index_references
            refname=sprintf('%s/%s%4.4i.edf',n1,refroot,n);

            % determine if reference needs to be created
            createref = ~checkreference(refname,sizenoheader);
            
            if createref  % reference does not exist or has wrong size and will be created
                if firstconversion
                    disp('Creating references ...')
                    firstconversion = 0;
                endif
                % check if single reference files exist
                if ~exist([n1 '/' refsname(nref,n,devserver)],'file')
                    n_used = refon*round(n/refon);
                else
                    n_used = n;
                endif

                pyhst_parfile(n1,'ref',nref,X,Y,1,'NO',100,0, \
                    'file_postfix',sprintf('_%4.4i.edf',n_used),'do_projection_median','YES','projection_median_filename',refpyhstname, \
                    'parfilename',refpyhstparfilename,'output_reconstruction','NO');
                if (exist ("OCTAVE_VERSION") == 5)
                    [status,output]=system(sprintf('pyhst %s/%s',n1,refpyhstparfilename));
                else
                    [output,status]=system(sprintf('pyhst %s/%s',n1,refpyhstparfilename));
                endif
                if status
                    error("Problems creating %s\n",refname);
                endif
                % store the reference file
	        ref = edfread(refpyhstname);
                edfwrite(refname,ref+0.5,datatype);
   	        printf('Reference %s%4.4i.edf written\r',refroot,n)

                % in case last reference is 999 we also create a link 1000 to this file
                if (n == nvue)
                    nround = refon*round(n/refon);
                    refnameround = sprintf('%s/%s%4.4i.edf',n1,refroot,nround);
                    if (n ~= nround) & (~exist(refnameround,'file'))
                        link(refname,refnameround);
                    endif
                endif
            endif % creation new ref file
        endfor
        unlink(refpyhstname);
        unlink(refpyhstparfilename);
        unlink('ref.vol');
    otherwise  % not save_ram
        f=zeros(X,Y,nref,'uint16');
%         f = cell(nref,1);
% attempt to read the data faster into a cell, but then data has to be reorganised for sorting
        ref=zeros(X,Y);
%         med = zeros(X,nref,'uint16');

        for n = index_references
            refname=sprintf('%s/%s%4.4i.edf',n1,refroot,n);

            % determine if reference needs to be created
            createref = ~checkreference(refname,sizenoheader);

            if createref
                if firstconversion
                    disp('Creating references ...')
                    firstconversion = 0;
                endif
                if isempty(find(index_refstobeconverted == n))
                    n_used = refon*round(n/refon);
                else
                    n_used = n;
                endif
    	        for k=1:nref
       	            name=[n1 '/' refsname(k,n_used,devserver)];
       	            f(:,:,k)=uint16(edfread(name));
    	        end
                whos('f')
    	        ref=median(f,3);
%     	        for k=1:nref
%        	            name=[n1 '/' refsname(k,n,devserver)];
%                     if ~exist(name,'file')
%                         name=[n1 '/' refsname(k,refon*round(n/refon),devserver)];
%                     endif
%        	            f{k}=uint16(edfread(name));
%     	        end
%                 whos('f')
%                 for l = 1:Y
%                     for k = 1:nref
%                         med(:,k) = f{k}(:,l);
%                     endfor
%                     %ref(:,l) = median(med,2);    
%                 endfor

	        edfwrite(refname,ref+0.5,datatype);
	        printf('Reference %s%4.4i.edf written\r',refroot,n)

                % in case last reference is 999 we also create a link 1000 to this file
                if (n == nvue)
                    nround = refon*round(n/refon);
                    refnameround = sprintf('%s/%s%4.4i.edf',n1,refroot,nround);
                    if (n ~= nround) & (~exist(refnameround,'file'))
                        link(refname,refnameround);
                    endif
                endif
            endif %createref
        endfor
    endswitch % (save_ram)
    clear('f')
else % average
    for n = index_references
        refname=sprintf('%s/%s%4.4i.edf',n1,refroot,n);

        % determine if reference needs to be created
        createref = ~checkreference(refname,sizenoheader);

        if createref
            if firstconversion
                disp('Creating references ...')
                firstconversion = 0;
            endif
            ref=zeros(X,Y);
            for k=1:nref
                name=[n1 '/' refsname(k,n,devserver)];
                if ~exist(name,'file')
                    name=[n1 '/' refsname(k,refon*round(n/refon),devserver)];
                endif
                ref += edfread(name);
            endfor

            ref /= nref;

            edfwrite(refname,ref+0.5,datatype);
            printf('Reference refHST%4.4i.edf written\r',n)

            % in case last reference is 999 we also create a link 1000 to this file
            if (n == nvue)
                nround = refon*round(n/refon);
                refnameround = sprintf('%s/%s%4.4i.edf',n1,refroot,nround);
                if (n ~= nround) & (~exist(refnameround,'file'))
                    link(refname,refnameround);
                endif
            endif
        endif %createref
    endfor
endif
if ~firstconversion
    printf('\n')
endif


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   Removing single reference files  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (deletesingleref)
    disp('Removing original reference files')
    for n = index_references
        refname=sprintf('%s/%s%4.4i.edf',n1,refroot,n);

        % determine if reference needs to be created
        createref = ~checkreference(refname,sizenoheader);

        if ~createref  % reference exists and has the good size
            if find(index_refstobeconverted == n)
	        % remove the single reference files
                for k=1:nref 
                    name=[n1 '/' refsname(k,n,devserver)];
	            unlink(name);
	        endfor
            endif
        else
	    printf('Reference %s%4.4i.edf is wrong, please check!!!\n',refroot,n)
        endif
    endfor
end

%------------------
function y=inputwdefault(question,default)
    lim=10;
    if length(default)>lim
        y=input(sprintf('%s      :\n[%s]\n',question,default),'s');
    else
        y=input(sprintf('%s      : [%s] ',question,default),'s');
    end
    if isempty(y)
        y=default;
    end
endfunction

%------------------
function s = checkreference(filename,sizenoheader)
    dirres = dir(filename);
    if isempty(dirres)
        s = 0;
    elseif (~dirres.bytes)
        s = 0;
    elseif (dirres.bytes ~= (length(header(filename))+sizenoheader))
        s = 0;
    else
        s = 1;
    endif
endfunction
