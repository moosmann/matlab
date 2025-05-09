function [free, available, total, cache] = free_memory()
% Host system memory
%
% [free, available, total, cache] = free_memory()

stats = [];
while isempty(stats)
    [~,w] = system('free -b | grep Mem');
    stats = str2double(regexp(w, '[0-9]*', 'match'));
    pause(0.5);
end
total = stats(1);
free = stats(3);
cache = stats(5);
available = stats(6);
