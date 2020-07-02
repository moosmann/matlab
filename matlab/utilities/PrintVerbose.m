function PrintVerbose(verbose, varargin)
% Print to standard output if 'verbose' is true, else do nothing.
%
% Written by Julian Moosmann.
% First version: 2016-09-28. Last modification: 2016-09-28
%
% function PrintVerbose(verbose, varargin)

%% Default
if nargin < 1
    verbose = 0;
end

%% Main
if verbose
    fprintf(1, varargin{:});
end