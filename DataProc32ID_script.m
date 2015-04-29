% Loop data reconstruction file of October data from 32ID@APS

%linfil = {[0 1],[1 1]};
%linfil = {[0 0],[1 0],[0 1],[1 1]}
%linfil = {[0 0 1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[1 1 1]};
linfil = {[1 1 0]};

for ff = 1:numel(linfil)
    DataProc32ID(1,1,linfil{ff},0);
end