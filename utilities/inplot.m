% function [h]=inplot(xvec,yvec,sran,c1,bds1,c2,bds2,c3,bds3)
%
% This function plots the data (xvec,yvec) and then plots a subset in the upper
% left corner of the screen to give you a "picture-in-a-picture".   
% Jeff Butera, jvbutera@eos.ncsu.edu, 30 November 95 (with much help from 
% Mike Buksas).
%
% Required INPUTS:
%   xvec - vector of x data points
%   yvec - vector of y data points, same length as xvec
%   sran - 2 vector of form [smin,smax] for subset to be plotted.
%          Note that (smin,smax) must be a subset of (min(xvec),max(xvec)).
%
% Optional INPUTS:
%   c   - a character, either 'w','p', or 's'
%   bds - 4 vector in form [xmin xmax ymin ymax]
%
%  If c=='p' then bds specifies the axis bounds for the outermost Plot.
%  If c=='s' then bds specifies the axis bounds for the Subplot.
%  If c=='w' then bds specifies the absolute bounds of the Window location.
%   The bounds for the window location will be "clipped" so the subplot window
%   does not run off the top or sides of the main plot.
%
% OUTPUTS:
%   h - a 2 vector of handles.  The first entry is the main graph and
%       the second is the subplot graph.
%
% Examples: x=[-500:500]; y=randn(1,1001);
%  inplot(x,y,[-200 -180])
%  inplot(x,y,[200 230],'w',[-400 400 4 10])
%  inplot(x,y,[300 315],'p',[-500 500 -3 8],'w',[-400 400 4 7])
%  inplot(x,y,[180 197],'s',[175 200 -1 1],'w',[-400 400 4 11])
function [h]=inplot(xvec,yvec,sran,c1,bds1,c2,bds2,c3,bds3)
%
ymax=max(yvec);
ymin=min(yvec);
xmax=max(xvec);
xmin=min(xvec);
%
% Figure out if user set window, entire plot or subplot axis, and which 
% is which...
%
wbds=[0 0 0 0];
pbds=[0 0 0 0];
sbds=[0 0 0 0];
if (nargin>3),
   if ((c1=='w')|(c1=='W')), 
     wbds=bds1;
   elseif ((c1=='p')|(c1=='P')),
     pbds=bds1;
   elseif ((c1=='s')|(c1=='S')),
     sbds=bds1;
   end
end
if (nargin>5),
   if ((c2=='w')|(c2=='W')), 
     wbds=bds2;
   elseif ((c2=='p')|(c2=='P')),
     pbds=bds2;
   elseif ((c2=='s')|(c2=='S')),
     sbds=bds2;
   end
end
if (nargin>7),
   if ((c3=='w')|(c3=='W')), 
     wbds=bds3;
   elseif ((c3=='p')|(c3=='P')),
     pbds=bds3;
   elseif ((c3=='s')|(c3=='S')),
     sbds=bds3;
   end
end
if (pbds==[0 0 0 0]),
   pbds=[xmin xmax ymin-0.15*(ymax-ymin) ymax+1.25*(ymax-ymin)];
end
plenx=pbds(2)-pbds(1);
pleny=pbds(4)-pbds(3);
%
% plot data in large window with plot axis already set
%
clf
plot(xvec,yvec)
axis(pbds)
hold on
%
% Check to see if user set window size manually
%
if (wbds==[0 0 0 0]),
   wbds=[0.2 0.6 0.3 0.25];
else
%
% Clip subplot window so it doesn't run off top or sides of figure
%
   tbds=wbds;
   wbds(1)=max(0.775*(tbds(1)-pbds(1))/plenx+0.13,0.2);
   wbds(2)=0.815*(max(tbds(3),pbds(1))-pbds(3))/pleny+0.11;  
   wbds(3)=min((0.775*(tbds(2)-pbds(1))/plenx)+0.13-wbds(1),0.85-wbds(1)); 
   wbds(4)=min((0.815*(tbds(4)-pbds(3))/pleny)+0.11-wbds(2),0.85-wbds(2));
end
%
% See if limits are valid
%
if (sran(1)>=sran(2)),
   disp(' ')
   disp(['Lower bound ',num2str(sran(1)),' must be less than upper bound ',num2str(sran(2))]);
   disp(' ')
   sran(1)=sran(2)-1;
end
%
% Clip data
%
i=min(find(xvec>=sran(1)));
j=max(find(xvec<=sran(2)));
%
% Check to see if user set window bounds manually...
%
if (sbds==[0 0 0 0]),
  sbds=[sran(1) sran(2) min(yvec(i:j)) max(yvec(i:j))];
else
  sbds(1)=min(sran(1),sbds(1));
  sbds(2)=max(sran(2),sbds(2));
end
%
% Plot delimiters...
%
sl=pbds(1)+plenx*(wbds(1)-0.13+wbds(3)*(sran(1)-sbds(1))/(sbds(2)-sbds(1)))/0.775;
sr=pbds(1)+plenx*(wbds(1)-0.13+wbds(3)*(sran(2)-sbds(1))/(sbds(2)-sbds(1)))/0.775;
y1=min(ymax,pbds(3)+pleny*((wbds(2)-0.11)/0.815-0.1));
y2=pbds(3)+pleny*((wbds(2)-0.11)/0.815-0.07);
y3=pbds(3)+pleny*((wbds(2)-0.11)/0.815-0.05);
plot([sran(1) sran(1) sl sl],[pbds(3) y1 y2 y3],'r-')
plot([sran(2) sran(2) sr sr],[pbds(3) y1 y2 y3],'r-')
hold off
%
% Switch axes handle for small plot...
%
gca1=gca;
gca2=axes('position',wbds);
plot(xvec,yvec)
axis(sbds)
axes(gca1);
if (nargout>0),
  h=[gca1 gca2];
end
%
return

