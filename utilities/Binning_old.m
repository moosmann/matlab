function [matrix,tbin]=Binning_old(matrix,resolution_binned,print_parameters)
    
if (nargin<3),print_parameters=0;end;

% Parameters.
[resx_unbinned,resy_unbinned] = size(matrix);
resx_binned                   = resolution_binned(1);
resy_binned                   = resolution_binned(2);
len                           = resx_unbinned*resy_unbinned;
binx                          = resx_unbinned/resx_binned; 
biny                          = resy_unbinned/resy_binned; 

% Print parameters.
if print_parameters==1,
fprintf('[resolution_unbinned]=[%g,%g], [resolution_binned]=[%g,%g], len=%g, [binx,biny]=[%g,%g]\n', ...
        resx_unbinned,resy_unbinned,resx_binned,resy_binned,len,binx,biny);
end;

% Binning.
tic;
if (binx==1)&&(biny==1),
    if print_parameters, fprintf(1,'No Binning\n'); end;
elseif (binx==2)&&(biny==2),
    matrix = (matrix(1:2:end,1:2:end)+matrix(2:2:end,1:2:end)+ ...
        matrix(1:2:end,2:2:end)+matrix(2:2:end,2:2:end))/4;
elseif (binx==4)&&(biny==4),
    matrix = (matrix(1:2:end,1:2:end)+matrix(2:2:end,1:2:end)+ ...
        matrix(1:2:end,2:2:end)+matrix(2:2:end,2:2:end))/4;
    matrix = (matrix(1:2:end,1:2:end)+matrix(2:2:end,1:2:end)+ ...
        matrix(1:2:end,2:2:end)+matrix(2:2:end,2:2:end))/4;
else
    matrix = reshape(sum(reshape(reshape(sum(reshape(matrix,binx,len/binx)),resx_binned,resy_unbinned)',biny,len/(binx*biny))),resx_binned,resy_binned)'/(binx*biny);
    if print_parameters, fprintf(1,'Binning method: reshape and sum\n'); end;
end;
tbin = toc;

% Flexible binning.
%intensity_down = zeros(resx_unbinned/xover,resy_unbinned/yover);
%for xx=1:xover, for yy=1:yover,  intensity_down = intensity_down +
%intensity(xx:xover:end,yy:yover:end); end; end;
%intensity = intensity_down/xover/yover; clear intensity_down;

