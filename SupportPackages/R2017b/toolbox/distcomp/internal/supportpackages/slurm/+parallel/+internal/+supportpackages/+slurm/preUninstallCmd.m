function preUninstallCmd()
%PREUNINSTALLCMD - internal: pre uninstall hook for Slurm

% Copyright 2015-2016 The MathWorks, Inc.
parallel.internal.supportpackages.preUninstallCmd('Slurm', 'parallel/slurm');
end
