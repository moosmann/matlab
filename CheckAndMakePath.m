function CheckAndMakePath(StringOfPath)
% Check if path 'StringOfPath' exists and if not create it.

if ~exist(StringOfPath,'dir')
    mkdir(StringOfPath);
end