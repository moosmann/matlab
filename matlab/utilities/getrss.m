function [m, p] = getrss( p )
% Returns memory usage of current MATLAB process as resident set size, the
% non-swapped physical memory that a task has used (in kiloBytes).


if nargin < 1
    p = feature('getpid');
end

%s = sprintf( 'ps -p %u u', p);
s = sprintf( 'ps -p %u -o rss --no-headers', p);
[~,r] = system(s);

m = str2double(r) / 1024;