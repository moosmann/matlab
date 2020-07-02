function figplot(varargin)
% Create figure with name optionally given by varargin{1} and forward input
% arguments to plot function.
%
% Written by Julian Moosmann, 2014-03-26, last version: 2014-03-26
%
% figplot(varargin)

if ischar(varargin{1})
    figure('Name',varargin{1})
    plot(varargin{2:end})
else
    figure
    plot(varargin{:})
end