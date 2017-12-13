function out = callerIsBaseWorkspace()

out = length( dbstack) == 1;
fprintf('\n caller is base workspace : %u\n', out )
