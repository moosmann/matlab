function [p,f] = last_folder_modified( scan_path )
% Returns path of the folder that was modified last
%
% ARGUMENTS
% scan_path: string, parent path or path pattern using wildcards.
%
% RETURNS
% p: string, absolute path to last folder modified
% f: string, folder last modified


if  nargin < 1
    scan_path = '/asap3/petra3/gpfs/p07/2021/data/11011705/raw/*';
end

d = dir( scan_path );
% fprintf('\n%u elements found', numel(d) )
ind = ones([numel(d),1],'logical');
for n = 1:numel(d)
    s = [d(n).folder filesep d(n).name];
    if strcmp( d(n).name(end), '.')        
        ind(n) = 0;
    end
    if ~isfolder( s )
       ind(n) = 0;
    end
end
d = d(ind);
% fprintf('\n%u folders found', numel(d) )
[~,di] = sort( [d(:).datenum] );
n = di(end);
f = d(n).name;
p = [d(n).folder filesep d(n).name];
