function varsize(Variable,Unit)
% Print memory allocation in kB, MB, or GB.

if nargin < 2
    Unit = 'MB';
end

switch lower(Unit)
    case {'kb','kilo','kilobyte'}
        mm = 1;
        Unit = 'kB';
    case {'mb','mega','megabyte'}
        mm = 2;
        Unit = 'MB';
    case {'gb','giga','gigabyte'}
        mm = 3;
        Unit = 'GB';
end

m = whos('Variable');
fprintf('\n %s   %g %s\n',inputname(1),m.bytes*1024^-mm,Unit)