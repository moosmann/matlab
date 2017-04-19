function [free, available, total, cache] = free_memory()
% Host system memory
%
% [free, available, total, cache] = free_memory()

[~,w] = unix('free -b | grep Mem');
stats = str2double(regexp(w, '[0-9]*', 'match'));
total = stats(1);
free = stats(3);
cache = stats(5);
available = stats(6);
