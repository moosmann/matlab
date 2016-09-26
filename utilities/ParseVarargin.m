function ParseVarargin(varargin)
% Assign input arguments defined in varargin in the workspace of the
% function that calls 'ParseVarargin'.
%
% Written by Julian Moosmann, last modified: 2013-09-06

varargin = varargin{:};
for nn = 1:2:numel(varargin)
    assignin('caller',varargin{nn},varargin{nn+1});
end