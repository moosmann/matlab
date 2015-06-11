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

## findshift
## function c = findshift(z1,z2,varargin)
##      returns the shift (2D) in pixels between two images
##      looks for a local optimum around the initial shift (see argument 4)
##      looks within a window given by argument 3
##      adapts the shift estimate in case optimum is at the edge of the window
##      result is c=[row column]
##      positive values corresponds to moving z2 to higher values of the index
##      to compensate drift: interpolate(f)(z2,row,column)
##      arguments:
##      argument 1: source image
##      argument 2: target image
##      argument 3 (optional): half window size to search for optimum shift(default: [5 5])
##      argument 4 (optional): approximative shift (default: [0 0])
##                  can be  [#rows #columns]
##                  or      {'fft','auto','fourier'} -> Fourier correlation
##                  or      {'manual','man','m'} -> manual (imagej_align)
##      argument 5 (optional): refine the result (default: 1)
##      e.g.    findshift(z1,z2)
##              findshift(z1,z2,2) -> searches in a window +/- 2
##              findshift(z1,z2,[10 2]) -> searches in a window +/- 10 in first direction, +/- 2 in second direction
##              findshift(z1,z2,1,'auto') -> takes as initial shift the result of correlate, uses a small window +/- 1
##              findshift(z1,z2,[],'man') -> take as initial shift the result of imagej_align
##              findshift(z1,z2,[],[],0) -> no refinement of the result
##              findshift(z1,z2,0,'man',0) -> result of manual alignment as it is
##              findshift(z1,z2,0,'fft',0) -> result of correlate as it is
##      see also correlate, imagej_align

## Author: P. Cloetens cloetens@esrf.fr
## 
## 2007-01-05 P. Cloetens cloetens@esrf.fr
## * Initial revision
## * based on correlaterealspace
## 2007-10-19 PC
## * Bug in help found by WL
## * third argument should be 0 or [0 0] to avoid real space correlation

function c = findshift(z1,z2,varargin)


if nargin < 2
    help findshift
    return
endif

if (nargin >= 3)
    if (length(varargin{1}) == 2)
        maxsize = varargin{1};
    else
        maxsize = [varargin{1} varargin{1}];
    endif
else
    maxsize = [];
endif
if isempty(maxsize)
    maxsize = [5 5];
endif

if (nargin >= 4)
    rapp = varargin{2};
else
    rapp = [];
endif
if isempty(rapp)
    rapp = [0 0];
endif

if (nargin >= 5)
    refine = varargin{3};
else
    refine = [];
endif
if isempty(refine)
    refine = 1;
endif


####################################
# determination of approximative shift (manually or Fourier correlation)

if !isnumeric(rapp)
    switch rapp
    case {'fft','auto','fourier'}
        rapp = correlate(z1,z2);
    case {'manual','man','m'}
        rapp = imagej_align(z1,z2);
    otherwise
        # user provides his own function to find approximative shift
        rapp = feval(rapp,z1,z2);
    endswitch
endif



####################################
# check if refinement with realspace correlation is required
# otherwise keep result as it is

if isequal(maxsize,[0 0])
    shiftfound = 1;
    if refine
        c = round(rapp);
    else
        c = rapp;
    endif
else
    shiftfound = 0;
    rapp = round(rapp);
endif

while (!shiftfound)
    maxsize = min(maxsize,size(z1)-abs(rapp));
    
    z1beg = max(1+rapp+maxsize,[1 1]);
    z1end = size(z1)+min(rapp-maxsize,[0 0]);

    cc=zeros(2*maxsize+1);
    
    z2beg = [1 1] + (rapp+maxsize>0).*size(cc)+(rapp+maxsize<=0).*(size(z1)-z1end+z1beg);
    z2end = z2beg + z1end - z1beg;

    z1p=z1(z1beg(1):z1end(1),z1beg(2):z1end(2))(:);

    for k = 1:size(cc,1)
        for l = 1:size(cc,2)
            z2p=z1p-z2(z2beg(1)-k:z2end(1)-k,z2beg(2)-l:z2end(2)-l)(:);
            cc(k,l) = sumsq(z2p-mean(z2p));
        endfor
    endfor
    
    [a,b]=min(cc);
    [d,c(2)]=min(a);
    c(1)=b(c(2));

    if !sum((c == 1)|(c == size(cc)))
        # check that we are not at the edge of the region that was sampled
        c = c + min2par(cc(c(1)-1:c(1)+1,c(2)-1:c(2)+1));
        shiftfound = 1;
    endif
    c = c + z1beg - z2beg;
    if !shiftfound
        rapp = c;
        printf('Changing shift estimate: %4d %4d\r',rapp(1),rapp(2))
        maxsize = min(maxsize,size(z1)-abs(rapp));
        if sum(maxsize == 0)
            printf('\nEdge of image reached\n')
            refine = 0;
            shiftfound = 1;
        endif
    else
        printf('\n')
    endif    
endwhile

####################################
# refine result; useful when shifts are not integer values

while (refine--)
    disp('Refining solution ...')
    z2n = interpolate(z2,c(1),c(2));
    z1p = circshift(z1,-(c>0).*ceil(abs(c)))(2:end-ceil(abs(c(1)))-1,2:end-ceil(abs(c(2)))-1)(:);
    z2n = circshift(z2n,-(c>0).*ceil(abs(c)))(1:end-ceil(abs(c(1))),1:end-ceil(abs(c(2))));
    ccrefine = zeros(3);
    for k = 1:3
        for l = 1:3
            z2p=z1p-z2n(4-k:end-k+1,4-l:end-l+1)(:);
            ccrefine(k,l) = sumsq(z2p-mean(z2p));
        endfor
    endfor
    crefine = min2par(ccrefine);    
    if abs(crefine) < 1
        c = c + crefine;
    else
        printf('Problems refining result\n')
    endif
endwhile


# subfunctions
#-------------------------
function mi=min2par(arg)

[X,Y]=meshgrid(-1:1,-1:1);
r=Y(:);
c=X(:);
xt=quad2func(r,c);
a=xt\arg(:);
mi=[-a(2) -a(3)]/[2*a(5) a(4);a(4) 2*a(6)];

#-------------------------
# function xt=quad2func(x,y)
# quadratic function of 2 parameters x and y
# size(xt)=[n 6] with n the number of observations
# for functions linear in the parameters:
# z = xt*a with size(a)=[6 1]
# i.e. z=a(1)+a(2)*x+a(3)*y+a(4)*x.*y+a(5)*x.^2+a(6)*y.^2
# if one has n observations given by z (size(z)=[n 1])
# the least squares estimation of a is
# a=xt\z;

function xt=quad2func(x,y)

xt=[ones(size(x)) x y x.*y x.^2 y.^2];
