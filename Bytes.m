function out = Bytes(array,PowerOf1024)
% Output (number of bytes allocated by input array.

if nargin < 2
    PowerOf1024 = 0;
end

out = whos('array');
out = out.bytes/1024^PowerOf1024;
