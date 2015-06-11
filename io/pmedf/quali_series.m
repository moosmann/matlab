## Copyright (C) 2006 P. Cloetens
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

## quali_series
## function retval = quali_series(direc,varargin)
##      checks for the motion during a scan
##      corrtot contains correlation between image before and after scan at angle 0
##      corrtot90 contains correlation between image at angle 90 during and after scan (total error can be double)
##      results are stored in quali.mat file
##      corr_imagesafterscan contains all correlations
##          (index 1: index in rot_positions of angle, index 2: horizontal/vertical)
##          the first dimension corresponds to the horizontal direction
##          the second dimension corresponds to the vertical direction
##      rot_projections contains angles (multiples of 90 degrees) for which correlation during/after scan is determined
##      rot_projections_index contains indices of projections corresponding to angles rot_projections
##      can process scans before or after flatfield correction
##      scanrange can be 180 or 360 degrees
##      images after scan can be acquired with or without inversion of rotation direction
##      results probably not very accurate when Scan_Type is continuous
##
##      arguments:
##      argument 1 (optional): directory name (absolute or relative) with scan
##                              ( default: current directory ) 
##      argument 2 (optional): method for shift determination ( default 5 or global FT_CORRELATEMETHOD )
##          0: only fft method 
##          n > 0: realspace correlation (findshift) of half window size correlatemethod 
##          n < 0: manual alignment (raw result)
##      argument 3 (optional): restrict to central part for 'accurate' method ( default 1 or global FT_CENTRALPART )
##
##      examples:
##          quali_series            -> processes current directory with real space correlation
##          quali_series(myscan)    -> processes myscan with real space correlation
##          quali_series([],-1)     -> processes current directory manually
##          quali_series([],0)      -> processes current directory with fft method
##          quali_series([],[],0)   -> processes current directory with 'accurate' method on complete image (slow)

## Author: P. Cloetens <cloetens@esrf.fr>
##
## 2006-08-02 P. Cloetens <cloetens@esrf.fr>
## * Initial revision
## * accepts scans that are flatfield corrected, choice of correlation method
## * if no *.info exists -> takes the last image number and subtracts 2
## 2007-02-23 P. Cloetens <cloetens@esrf.fr>
## * use findshift instead of correlate + correlaterealspace
## * add size in printout (to determine plane with minimum shift in holotomo_slave)
## 2008-01-14 P. Cloetens <cloetens@esrf.fr>
## * rewrite to adapt to scan over 360 degrees
## * possibility to determine the shifts manually (correlatemethod negative)
## * stores corrtot and corrtot90 as before, but also corr_imagesafterscan (different angles)
## 2008-01-15 P. Cloetens <cloetens@esrf.fr>
## * determination inversion/no inversion at end of scan improved
## 2008-01-18 P. Cloetens <cloetens@esrf.fr>
## * do not restrict to central part in case of manual method
## * add third argument centralpart and global variable FT_CENTRALPART
## * add warning when Scan_Type is continuous and ANGLE_BETWEEN_PROJECTIONS has been adjusted manually
## * asks for nvue when info-file does not exist
## 2008-01-18 PC
## * remove bug that prevented to run quali_series outside the directory to be processed
## * remove create_projection as sub-function, make it an independent function

function varargout = quali_series(direc,correlatemethod,centralpart)
    global FT_CORRELATEMETHOD FT_CENTRALPART
 
    mot_names_rot_motor = {'pmo','rots','rl'}; # list with names of tomography motor on different instruments

    if !exist('correlatemethod','var')
        correlatemethod = [];
    endif
    if !exist('centralpart','var')
        centralpart = [];
    endif
    if !exist('direc','var')
        direc = '';
    endif
    if isempty(direc)
        direc = pwd;       
    endif
    
    if isempty(correlatemethod)
        if isempty(FT_CORRELATEMETHOD)
            correlatemethod = 5;
        else
            correlatemethod = FT_CORRELATEMETHOD;
        endif
    endif
    if isempty(centralpart)
        if isempty(FT_CENTRALPART)
            centralpart = 1;
        else
            centralpart = FT_CENTRALPART;
        endif
    endif

    disp(['DIRECTORY IS ' direc])

    ################################################################
    ##################   Information on job   ######################
    ################################################################

    # interactive information
    n1 = cleandirectoryname(direc);

    # n2: fileprefix, taken from directory name
    pos = findstr(n1,'/');
    if isempty(pos)
    	n2 = n1;
    elseif	(pos(end) == length(n1))
    	n2 = n1(pos(end-1)+1:end-1);
    else
    	n2 = n1(pos(end)+1:end);
    endif
    	
    # information from edf-files
    fname = [n1 '/' n2 '0000.edf'];
    if !exist(fname,"file")
    	n2 = input('name of projections   : ','s');
    	fname = [n1 '/' n2 '0000.edf'];
    endif
    hd = header(fname);
    X = findheader(hd,'Dim_1','integer');
    Y = findheader(hd,'Dim_2','integer');
    pixsize = findheader(hd,'optic_used','float');

    # information from info-file
    infofilename = [n1 '/' n2,'.info'];
    if exist(infofilename,"file")  # *.info exists
        fp = fopen(infofilename);
        hd = fscanf(fp,'%c');
        fclose(fp);
        if isempty(pixsize)
            pixsize = findheader(hd,'Optic＿used','float');
        endif
        nvue = findheader(hd,'TOMO_N','integer');
        refon = findheader(hd,'REF_ON','integer');
        nref = findheader(hd,'REF_N','integer');
        ndark = findheader(hd,'DARK_N','integer');
        scanrange = findheader(hd,'ScanRange','float');
        Scan_Type = findheader(hd,'Scan_Type','string');
    else	# *.info does not exist
        # we assume that flatfield is already done and *.info was not copied/created
        pixsize = [];
        nvue = input('Number of projections   : ');
        refon = [];
        scanrange = [];
    endif
    if isempty(nvue)
        # we take the last image number and subtract 2
        file_list = glob(sprintf('%s/%s????.edf',n1,n2));
        nvue = str2num(file_list{end}(length(n1)+length(n2)+2:end-4))-scanrange/90;
    endif
    if isempty(scanrange)
        # check if scanrange was 180 or 360 degrees
        # not very nice, but maybe info was not in info-file or xml-file
        if exist(sprintf('%s/%s%4.4i.edf',n1,n2,nvue+4),'file')
            is360 = 1;
            scanrange = 360;
        else
            is360 = 0;
            scanrange = 180;
        endif
    else
        is360 = (scanrange == 360);
    endif
    parfilename = [n1 '/' n2 'slice.par'];
    if exist(parfilename,"file")
        # read angle_between_projections
        fid = fopen(parfilename);
        hd = fscanf(fid,'%c');
        fclose(fid);
        angle_between_projections = findheader(hd,'ANGLE_BETWEEN_PROJECTIONS','float32');
    else
        angle_between_projections = scanrange/nvue;
    endif
    # display some conclusions and warnings
    printf('ScanRange was %d degrees\n',scanrange)
    if (abs(angle_between_projections - scanrange/nvue) > 1e-5)
        printf("ANGLE_BETWEEN_PROJECTIONS = %0.5f according to slice.par file\n",angle_between_projections)
        printf("ANGLE_BETWEEN_PROJECTIONS = %0.5f according to number of projections\n",scanrange/nvue)
        disp("!!! QUALI_SERIES MAY GIVE POOR RESULTS !!!")
    endif
    if !isempty(strmatch('CONT',upper(Scan_Type)))
        printf("Scan_Type was %s\n",Scan_Type)
        disp("!!! QUALI_SERIES MAY GIVE POOR RESULTS !!!")
    endif
    # determine if images after scan are taken in opposite rotation direction (case ID19) or same (case ID22NI)
    rot_pos_imageafterscan = [];
    k = 1;
    while ( isempty(rot_pos_imageafterscan) & (k <= length(mot_names_rot_motor)) )
        rot_pos_imageafterscan = header(sprintf('%s/%s%4.4i.edf',n1,n2,nvue+1),mot_names_rot_motor{k++},'motor');
    endwhile
    rot_positions = 0:90:scanrange-90;
    rot_pos_imageafterscan = mod(rot_pos_imageafterscan, 360); # in case motor position is not modulo 360
    if isempty(rot_pos_imageafterscan)
        disp("Motor position can not be read")
        rot_positions_index = (scanrange/90:-1:1);
    elseif ( abs(rot_pos_imageafterscan - (scanrange-90)) < 0.5 )
        rot_positions_index = (scanrange/90:-1:1); # inversion
    elseif ( is360 & ( abs(rot_pos_imageafterscan -90) < 0.5 ) )
        rot_positions_index = [4 1:3]; # no inversion
    elseif ( !is360 & ( abs(header(sprintf('%s/%s%4.4i.edf',n1,n2,nvue+2),mot_names_rot_motor{k-1},'motor') -90) < 0.5 ) )
        rot_positions_index = [1 2]; # no inversion, we check for image nvue+2 because nvue+1 tells little
    else
        disp("Motor position is not consistent")
        rot_positions_index = (scanrange/90:-1:1);
    endif
    rot_positions_index += nvue;
    # display assumed angles as check
    printf("Assuming projection ")
    for k = 1:length(rot_positions_index)
        printf("%d deg: %d; ",rot_positions(k),rot_positions_index(k))
    endfor
    printf("\n")
    
    # test for refHST
    do_flatfield = 1;
    if exist(sprintf('%s/refHST%04i.edf',n1,0))
        refstring = 'refHST';
    elseif exist(sprintf('%s/ref0000_%04i.edf',n1,0))
        refstring = 'ref0000_';
    elseif exist(sprintf('%s/ref1%04i.edf',n1,0))
        refstring = 'ref1';
    else
        do_flatfield = 0;
        refstring = [];
    endif

    ################################################################
    ##################   Creating dark file   ######################
    ################################################################

    if do_flatfield
        ###########
        # read dark.edf or create if not existing
        fname = [n1 '/dark.edf'];
        if exist(fname,"file")
            d = edfread(fname);
        else
            fname_summed = [n1 '/darkend0000.edf'];
            if isempty(ndark)
                ndark = header(fname_summed,'no_summed','integer');
                if isempty(ndark)
            		disp('Number of dark images is not specified, exiting !!!')
            		return
                endif
            endif
    	    d = edfread(fname_summed)/ndark;
    	    edfwrite(fname,d+0.5,'uint16');
        endif
    else
        d = [];
    endif

    ################################################################
    #######   Setting ROI depending on correlation method   ########
    ################################################################
    
    if ( (correlatemethod > 0) & centralpart )
        cuth = min(1024,round(X/2));
        cutv = min(1024,round(Y/2));
    else
        cuth = X;
        cutv = Y;
    endif
    printf("Using image size: %d (h) x %d (v)\n",cuth,cutv)
    c2=[round((X-cuth)/2)+1 round((X-cuth)/2)+cuth round((Y-cutv)/2)+1 round((Y-cutv)/2)+cutv];

    ####################################################################################
    ##################   Create projections and determine shifts   #####################
    ####################################################################################
    
    corr_imagesafterscan = zeros(length(rot_positions),2);
    for k = 1:length(rot_positions)
        # projection during scan
        a0 = create_projection(n1,n2,nvue*(rot_positions(k)/scanrange),do_flatfield,refstring,refon,d);
        # projection after scan
        a0e = create_projection(n1,n2,rot_positions_index(k),do_flatfield,refstring,refon,d);
        # shift determination
        if (correlatemethod > 0)
            # real space correlation
            corr = findshift(a0(c2(1):c2(2),c2(3):c2(4)),a0e(c2(1):c2(2),c2(3):c2(4)),correlatemethod);
        elseif (correlatemethod < 0)
            # manual correlation
            corr = findshift(a0(c2(1):c2(2),c2(3):c2(4)),a0e(c2(1):c2(2),c2(3):c2(4)),0,'man',0);
        else
            # correlation using fft
            corr = correlate(a0(c2(1):c2(2),c2(3):c2(4)),a0e(c2(1):c2(2),c2(3):c2(4)));
        endif
        printf('Correlation at %d degrees during/after scan: %+0.4f (h) %+0.4f (v); size: %+0.4f\n', \
                rot_positions(k),corr(1),corr(2),sqrt(corr*corr'))
        corr_imagesafterscan(k,:) = corr;
        # find angles at 0 and 90 degrees for backward compatibility
        switch rot_positions(k)
            case 0
                corrtot = corr;
            case 90
                corrtot90 = corr;
        endswitch 
    endfor
    
    ####################################################################################
    ##########################   Save and output results   #############################
    ####################################################################################
    
    save([n1,'/quali.mat'],'corrtot','corrtot90','corr_imagesafterscan','rot_positions','rot_positions_index')

    switch nargout
        case 1
            varargout{1} = [corrtot;corrtot90];
        case 2
            varargout{1} = corrtot;
            varargout{2} = corrtot90;
    endswitch
endfunction
