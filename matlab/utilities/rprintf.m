function rprintf(varargin)
% Print and log function for reco script printing to reco log file and
% command lin.

global fid_reco
fprintf(fid_reco,varargin{:});
fprintf(varargin{:})
