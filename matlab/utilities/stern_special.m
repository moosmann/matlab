function c = stern_special
%

%   Copyright 1984-2015 The MathWorks, Inc.
% 
% if nargin < 1
%    f = get(groot,'CurrentFigure');
%    if isempty(f)
%       m = size(get(groot,'DefaultFigureColormap'),1);
%    else
%       m = size(f.Colormap,1);
%    end
% end

load('cm_stern_special.mat','cm_stern_special');
c = cm_stern_special;
