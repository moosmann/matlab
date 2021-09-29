function [m, p] = getmem( p )
% Returns memory usage in percent of current MATLAB process as ratio of the
% process's resident set size  to the physical memory on the machine,
% expressed as a percentage.  


if nargin < 1
    p = feature('getpid');
end

%s = sprintf( 'ps -p %u u', p);
s = sprintf( 'ps -p %u -o %%mem --no-headers', p);
[~,r] = system(s);

m = str2double(r);